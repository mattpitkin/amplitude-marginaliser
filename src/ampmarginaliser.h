/*
*  Copyright (C) 2014 Matthew Pitkin
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/** \file ampmarginaliser.h
 * \brief Amplitude Marginaliser functions
 * 
 * This library of functions will analytically marginalise a Gaussian likelihood ratio
 * (i.e. the ratio of a Gaussian likelihood given a model to a Gaussian likelihood
 * for pure noise) over the amplitudes (e.g., the five amplitudes \f$A\f$, \f$B\f$, \f$C\f$,
 * \f$D\f$ and \f$E\f$) of an arbitrary number of independent components of a signal model
 * (e.g., of the five model components \f$Af(x_1)\f$, \f$Bg(x_2)\f$, \f$Ch(x_3)\f$,
 * \f$Di(x_4)\f$ and \f$Ej(x_5)\f$):
   \f[
   \mathcal{L}(x_1, x_2, x_3, x_4, x_5) = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}
   \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \int_{-\infty, 0}^{\infty}
   \exp{-\sum_{n=1}^N\frac{1}{2\sigma^2} (Af(x_1)_n + Bg(x_2)_n + Ch(x_3)_n + Di(x_4)_n +
   Ej(x_5)_n)^2 - 2d_n(Af(x_1)_n + Bg(x_2)_n + Ch(x_3)_n + Di(x_4)_n + Ej(x_5)_n) ) }
   {\rm d}A {\rm d}B {\rm d}C {\rm d}D {\rm d}E.
   \f]
 * In the above equation we leave it open as to whether the final model amplitude is
 * marginalised between \f$-\infty\f$ and \f$\infty\f$, or 0 and \f$\infty\f$ (i.e. is only 
 * allowed to be positive). When expanded out and terms grouped such that summations on
 * pairs of models, and models on data, are done separately the equation for the first
 * integral e.g. over \f$A\f$ the equation can be placed into the format
 \f[
 \mathcal{L}_A = \int_{-\infty}^{\infty} \exp{\left(-\frac{1}{2\sigma^2}(XA^2 + YA + Z)\right)}
 {\rm d}A, 
 \f]
 * which has the solution (via e.g. <A HREF="http://www.wolfram.com/mathematica/">Mathematica</A>)
 \f[
 \mathcal{L}_A = \sqrt{\frac{2\pi}{X}}\exp{\left(-\frac{1}{2\sigma^2}(Z-Y^2/(4X))\right)}.
 \f]
 * For that case that the integral is between 0 and $\infty$ this is instead
 \f[
 \mathcal{L}_A = \sqrt{\frac{\pi}{2X}}\exp{\left(-\frac{1}{2\sigma^2}(Z-Y^2/(4X))\right)}
 {\rm erfc}{\left(\frac{Y}{2\sqrt{2\sigma^2X}}\right)}.
 \f]
 * Mathematica can in fact solve the full equation (for five amplitudes at least), presumably
 * by repeated iterations of the above formula, but when the number of models is greater
 * than 3 the output is so long that it is unusable. However, it is possible to iterate the
 * above equation numerically as demonstrated in these library functions.
 * 
 * As this is a likelihood ratio you can get back to the likelihood by multiplying by:
 \f[
 (2\pi\sigma^2)^{-N/2}\exp{\left(-\sum_{n=1}^N \frac{d_n}{2\sigma^2} \right)}.
 \f]
 * 
 * These functions do not calculate a Bayes factor as no model parameter priors are used.
 * The functions integrate the model amplitudes over an infinite range and they don't
 * explicitly include prior values on them. However, if one can assign a flat prior (or
 * wide Gaussian prior) on the amplitudes, and the likelihood falls of to zero well
 * within that prior range, then amplitude priors can be subsequently applied to the output.
 * Priors on any other (non-independent amplitude) model parameters can be applied as
 * required as this likelihood ratio must still be calculated over those parameter ranges.
 * However, if actually comparing Bayesian evidences in which the same models have been
 * used then using this likelihood ratio is fine as the priors would cancel out.
 *
 * It should also be noted that as in <A HREF="http://bayes.wustl.edu/glb/book.pdf">Bretthorst</A>
 * the marginalisation integral could be
 * simplified by converting the independent model component into a new set of
 * orthogonal models, which makes the model component cross terms become zero. In the
 * orthogonal model space the the likelihood model is much simpler, but does require
 * pre-calculating the orthogonal models. Another issue with this is that it requires
 * that all the model amplitudes are integrated between \f$-\infty\f$ and \f$\infty\f$.
 */

#include <math.h>
#include <string.h>

/* GSL headers */
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define AM_SQUARE(x) (x*x)          /**< The square of \c x */
#define AM_LN2PI (M_LN2 + M_LNPI)   /**< \f$ \ln{(2\pi)} \f$ using <A HREF="https://www.gnu.org/software/gsl/manual/html_node/Mathematical-Constants.html">GSL constants</a> */
#define AM_LNPI_2 (M_LNPI - M_LN2)  /**< \f$ \ln{(\pi/2)} \f$ */

/* define functions */

/** \brief Marginalise a Gaussian likelihood ratio over the amplitudes for each of a set of
  * model components 
  * 
  * Function to marginalise over the amplitudes for each of the model components, where there can
  * be an arbitrary number of models. In general all the amplitudes are marginalised between
  * \f$-\infty\f$ and \f$\infty\f$, although if specified the final model component amplitude will be
  * marginalised between 0 and \f$\infty\f$ i.e. it's a purely positive model.
  *
  * @param Nmodels - the number of components of the signal model with separate amplitudes
  * @param modelModel - a 2D array containing the sums of the model cross products from the likelihood
  * e.g. for three model components \f$f\f$, \f$g\f$, \f$h\f$, this should contain:
  * \f[
  * \left(
  * \begin{array}{c c c}
  * \sum f_i^2 & \sum f_i g_i & \sum f_i h_i \\
  * 0 & \sum g_i^2 & \sum g_i h_i \\
  * 0 & 0 & \sum h_i^2
  * \end{array}
  * \right)
  \f]
  * @param dataModel - a 1D array containing the sums of the data (\f$d\f$) with each model component from the
  * likelihood e.g. \f${\sum_i d_i f_i), \sum d_i g_i, \sum d_i h_i}\f$.
  * @param sigma - the noise standard deviation of the data.
  * @param lastHalfRange - if this is set to 1 the final model component will be marginalised between
  * 0 and infinity rather than the default of -infinity to infinity.
  * 
  * @return logL - the natural logarithm of the marginalised likelihood ratio.
  * 
  * Note: this function does not apply priors to the marginalised amplitudes, but flat prior
  * values could be applied afterwards provided the likelihoods have approached zero towards
  * the edges of the prior range.
  */  
double marginalise_amplitudes(int Nmodels,  
                              double **modelModel, 
                              double *dataModel, 
                              double sigma,
                              unsigned int lastHalfRange);

/** \brief The same as \c marginalise_amplitudes(), but for \c float inputs rather than \c double.
 */
float marginalise_amplitudes_f(int Nmodels,
                               float **modelModel,
                               float *dataModel,
                               float sigma,
                               unsigned int lastHalfRange);

/** \brief The same as \c marginalise_amplitudes(), but with the \c modelModel array made into
 * a 1D vector rather than a 2D array.
 */
double marginalise_amplitudes_linear(int Nmodels,
                                     double modelModel[],
                                     double dataModel[],
                                     double sigma,
                                     unsigned int lastHalfRange);

/** \brief Performs the likelihood ratio marginalisation explicitly written out for the model amplitude components. 
 * 
 * For three amplitudes the marginalised likelihood ratio is given by:
 \f[
 \mathcal{L} = (2\pi)^{3/2}\sigma \frac{e^{-\frac{ d_h^2(f_g^2-FG) -2d_h(d_g f_g f_h - d_f f_h G -d_g F g_h + d_f f_g
g_h) + d_g^2(f_h^2 - FH) + 2d_f d_g(f_g H - f_h g_h) + d_f^2(g_h^2 - GH)}{2\sigma^2 (FGH + 2f_g f_h g_h -G f_h^2 
- H f_g^2 - F g_h^2) }}}{\sqrt{F}\sqrt{\frac{(FG - f_g^2)}{F}}\sqrt{\frac{FGH + 2f_g f_h g_h -G f_h^2 - F g_h^2 -
H f_g^2}{FG - f_g^2}}}
 \f]
 * for an integral between \f${-\infty, \infty}\f$ for all amplitudes, and
 \f{eqnarray*}{
   \mathcal{L} =& \sqrt{2}\pi^{3/2}\sigma \frac{e^{-\frac{ d_h^2(f_g^2-FG) -2d_h(d_g f_g f_h - d_f f_h G -d_g F g_h +
d_f f_g
g_h) + d_g^2(f_h^2 - FH) + 2d_f d_g(f_g H - f_h g_h) + d_f^2(g_h^2 - GH)}{2\sigma^2 (FGH + 2f_g f_h g_h -G f_h^2 
- H f_g^2 - F g_h^2) }}}{\sqrt{F}\sqrt{\frac{(FG - f_g^2)}{F}}\sqrt{\frac{FGH + 2f_g f_h g_h -G f_h^2 - F g_h^2 -
H f_g^2}{FG - f_g^2}}} \times \\
 & \left( 1 + {\rm erf}\left( \frac{d_f f_g g_h + d_g f_g f_h + FG d_h - G d_f f_h - F d_g g_h - d_g
f_g^2}{\sqrt{2}\sigma(FG-f_g^2)\sqrt{\frac{FGH + 2f_g f_h g_h - F g_h^2 - G f_h^2 - h f_g^2}{FG-f_g^2}}} \right)
\right)
 \f}
 * for the final integral being between $[0, \infty]$. In these equations we have used the following
 * substitutions: \f$d_f = \sum_{i=1}^N d_i f_i\f$,
 * \f$d_g = \sum_{i=1}^N d_i g_i\f$, \f$d_h = \sum_{i=1}^N d_i h_i\f$, \f$f_g = \sum_{i=1}^N f_i g_i\f$
 * \f$f_h = \sum_{i=1}^N f_i h_i\f$, \f$g_h = \sum_{i=1}^N g_i h_i\f$, \f$F = \sum_{i=1}^N f_i^2\f$, 
 * \f$G = \sum_{i=1}^N g_i^2\f$, and \f$H = \sum_{i=1}^N H_i^2\f$ where \f$d\f$ is
 * the data and \f$f\f$, \f$g\f$ and \f$h\f$ are the three model components for which the
 * amplitudes have been marginalised over.  

 * This function explicitly writes out the marginalised likelihood ratio for a model consisting
 * of three components with the amplitudes marginalised over. This can be used to test the correctness
 * of the \c marginalise_amplitudes() function. */
double marginalise_three_amplitudes(double **modelModel,
                                    double *dataModel,
                                    double sigma,
                                    unsigned int lastHalfRange);

/** \brief The same as \c marginalise_three_amplitudes(), but for \c float inputs rather than \c double.
 */
float marginalise_three_amplitudes_f(float **modelModel,
                                     float *dataModel,
                                     float sigma,
                                     unsigned int lastHalfRange);
