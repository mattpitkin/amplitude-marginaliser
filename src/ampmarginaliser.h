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
 * \brief A set of functions for analytically marginalising a Gaussian likelihood ratio
 * over signal model component amplitudes.
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
  * e.g. for three model components (see \ref sec_desc) \f$f\f$, \f$g\f$, \f$h\f$, this should contain:
  * \f[
  * \left(
  * \begin{array}{c c c}
  * \sum_i f_i^2 & \sum_i f_i g_i & \sum_i f_i h_i \\
  * 0 & \sum_i g_i^2 & \sum_i g_i h_i \\
  * 0 & 0 & \sum_i h_i^2
  * \end{array}
  * \right)
  \f]
  * @param dataModel - a 1D array containing the sums of the product of the data, \f$d\f$, with each model component from the
  * likelihood e.g. \f${\sum_i d_i f_i, \sum_i d_i g_i, \sum_i d_i h_i}\f$.
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

/** \brief Marginalise over all amplitudes except for the final model component.
 *
 * This function allows for a model that contains components for which the amplitudes
 * are required to be marginalisation over, and also a component (which could contain multiple
 * components defined by various parameters itself) for which no marginalisation is
 * required. It follows the same route as for \c marginalise_amplitudes(), but stops
 * before the final marginalisation.
 *
 * All the amplitude marginalisations are performed over the \f${-\infty, \infty}\f$ range.
 */
double marginalise_amplitudes_except_final(int Nmodels,
                                           double **modelModel,
                                           double *dataModel,
                                           double sigma);

/** \brief The same as \c marginalise_amplitudes_except_final(), but for \c float inputs rather than \c double.
 */
float marginalise_amplitudes_except_final_f(int Nmodels,
                                            float **modelModel,
                                            float *dataModel,
                                            float sigma);

/** \brief The same as \c marginalise_amplitudes(), but with the \c modelModel array made into
 * a 1D vector rather than a 2D array.
 */
double marginalise_amplitudes_linear(int Nmodels,
                                     double modelModel[],
                                     double dataModel[],
                                     double sigma,
                                     unsigned int lastHalfRange);

/** \brief The same as \c marginalise_amplitudes_except_final(), but with the \c modelModel array made into
 * a 1D vector rather than a 2D array.
 */
double marginalise_amplitudes_except_final_linear(int Nmodels,
                                                  double modelModel[],
                                                  double dataModel[],
                                                  double sigma);

/** \brief Performs the likelihood ratio marginalisation explicitly written out for the model amplitude components.
 *
 * For three amplitudes the marginalised likelihood ratio is given by:
 \f[
 \mathcal{L} = (2\pi)^{3/2}\sigma^3 \frac{e^{-\frac{ d_h^2(f_g^2-FG) -2d_h(d_g f_g f_h - d_f f_h G -d_g F g_h + d_f f_g
g_h) + d_g^2(f_h^2 - FH) + 2d_f d_g(f_g H - f_h g_h) + d_f^2(g_h^2 - GH)}{2\sigma^2 (FGH + 2f_g f_h g_h -G f_h^2
- H f_g^2 - F g_h^2) }}}{\sqrt{F}\sqrt{\frac{(FG - f_g^2)}{F}}\sqrt{\frac{FGH + 2f_g f_h g_h -G f_h^2 - F g_h^2 -
H f_g^2}{FG - f_g^2}}}
 \f]
 * for an integral between \f${-\infty, \infty}\f$ for all amplitudes, and
 \f{eqnarray*}{
   \mathcal{L} =& \sqrt{2}\pi^{3/2}\sigma^3 \frac{e^{-\frac{ d_h^2(f_g^2-FG) -2d_h(d_g f_g f_h - d_f f_h G -d_g F g_h +
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

/** \brief Explicitly marginalise over the first three model component amplitudes, but not the final one.
 *
 * This function uses the fully expanded integral over the amplitudes of first three components of a
 * signal model, whilst not integrating over the final fourth component's amplitude. All the amplitude
 * marginalisations are performed over the \f${-\infty, \infty}\f$ range. This can be used to test the
 * correctness of the \c marginalise_amplitudes_except_final() function.
 */
double marginalise_three_amplitudes_exclude_final(double **modelModel,
                                                  double *dataModel,
                                                  double sigma);

