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

#include "ampmarginaliser.h"

/** Function to marginalise over the amplitudes for each of the model components, where there can
  * be an arbitrary number of models. In general all the amplitudes are marginalised between
  * -infinity and infinity, although if specified the final model component amplitude will be
  * marginalised between 0 and infinity i.e. it's a purely positive model.
  *
  * Inputs:
  * \c Nmodels - the number of components of the signal model with separate amplitudes
  * \c modelModel - a 2D array containing the sums of the model cross products from the likelihood
  * e.g. for three model components f, g, h, this should contain:
  * [[ sum(f*f), sum(f*g), sum(f*h) ],
  *  [       0., sum(g*g), sum(g*h) ],
  *  [       0.,       0., sum(h*h) ]].
  * \c dataModel - a 1D array containing the sums of the data (d) with each model component from the
  * likelihood e.g. [ sum(d*f), sum(d*g), sum(d*h) ].
  * \c sigma - the noise standard deviation of the data.
  * \c lastHalfRange - if this is set to 1 the final model component will be marginalised between
  * 0 and infinity rather than the default of -infinity to infinity.
  * 
  * Returns:
  * \c logL - the natural logarithm of the marginalised likelihood.
  * 
  * Note: this function does not apply priors to the marginalised amplitudes, but flat prior
  * values could be applied afterwards provided the likelihoods have approached zero towards
  * the edges of the prior range.
  */  
double marginalise_amplitudes(int Nmodels, double **modelModel, double *dataModel, double sigma,
                              unsigned int lastHalfRange){
  /* coefficients of squares of model amplitudes */
  double squared[Nmodels];

  /* coefficients of model amplitudes */
  double coeffs[Nmodels][Nmodels];
  memset(coeffs, 0, sizeof(coeffs[0][0]) * Nmodels * Nmodels); /* initialise all values to zero */

  int i = 0, j = 0, k = 0, nm = Nmodels-1;

  double X = 0., invX = 0., invTwoX = 0., invFourX = 0., Y = 0., Z = 0., logL = 0.;

  /* set up coeffs matrix */
  for ( i=0; i<Nmodels; i++ ){
    squared[i] = modelModel[i][i];
    coeffs[i][i] = -2.*dataModel[i];

    if( isinf(squared[i]) || isinf(coeffs[i][i]) ){ return -INFINITY; }

    for ( j=(i+1); j<Nmodels; j++ ){
      coeffs[i][j] = 2.*modelModel[i][j];
     
      if( isinf(coeffs[i][j]) ){ return -INFINITY; }
    }
  }
  
  /* for five models the above coeff matrix would be equivalent to:
   * 
   * {{ -2.*DF,  2.*fg,  2.*fh,  2.*fi,  2.*fj},
   *  {     0., -2.*DG,  2.*gh,  2.*gi,  2.*gj},
   *  {     0.,     0., -2.*DH,  2.*hi,  2.*hj},
   *  {     0.,     0.,     0., -2.*DI,  2.*ij},
   *  {     0.,     0.,     0.,     0., -2.*DJ}}
   * where e.g. DF = sum(d*f) and fg = sum(f*g), and where for models with
   * amplitudes A, B, C, D and E these are the coefficients from the likelihood
   * for each, with rows corresponding to coefficients of each amplitude A
   * through E, and columns corresponding to cross term amplitudes, e.g.
   * coeffs[0][0] is terms with coefficients just of A,
   * coeffs[0][1] is terms with coefficients of A and B
   * coeffs[0][1] is terms with coefficients of A and C
   * ...
   * coeffs[2][2] is terms with coefficients of just C
   * coeffs[2][3] is terms with coefficients of C and D
   * ...
   * coeffs[4][4] is terms with coefficients just of E
   */

  for ( i=0; i<nm; i++ ){
    X = squared[i];
    invX = 1./X;
    invTwoX = 0.5*invX;
    invFourX = 0.25*invX;
    
    /* get the coefficients from the Y^2 term */
    for ( j=i; j<Nmodels; j++ ){
      for ( k=i; k<Nmodels; k++ ){
        /* add on new coefficients of squared terms */
        if ( j == i ){
          if ( k > j ){ squared[k] -= coeffs[j][k]*coeffs[j][k]*invFourX; }
          else if ( k == j ){ Z -= coeffs[j][k]*coeffs[j][k]*invFourX; }
        }
        else {
          if ( k == i ){ coeffs[j][j] -= coeffs[i][j]*coeffs[i][k]*invTwoX; }
          else if ( k > j ){ coeffs[j][k] -= coeffs[i][j]*coeffs[i][k]*invTwoX; }
        }
      }
    }
  }

  X = squared[nm];
  Y = coeffs[nm][nm];

  /* calculate analytic integral and get log likelihood */
  for ( i=0; i<Nmodels; i++ ){ logL -= 0.5*log(squared[i]); }

  logL -= 0.5*(Z - 0.25*Y*Y/X) / (sigma*sigma);

  logL += log(sigma);
  
  /* check whether final model is between 0 and infinity or -infinity and infinity */
  if (lastHalfRange == 1 ){
    logL += 0.5*(double)nm*AM_LN2PI + 0.5*AM_LNPI_2 + gsl_sf_log_erfc(0.5*Y/(sigma*sqrt(2.*X)));
  }
  else{ logL += 0.5*(double)Nmodels*AM_LN2PI; }

  return logL;
}


float marginalise_amplitudes_f(int Nmodels, float **modelModel, float *dataModel, float sigma,
                               unsigned int lastHalfRange){
  /* coefficients of squares of model amplitudes */
  float squared[Nmodels];

  /* coefficients of model amplitudes */
  float coeffs[Nmodels][Nmodels];
  memset(coeffs, 0, sizeof(coeffs[0][0]) * Nmodels * Nmodels); /* initialise all values to zero */

  int i = 0, j = 0, k = 0, nm = Nmodels-1;

  float X = 0., invX = 0., invTwoX = 0., invFourX = 0., Y = 0., Z = 0., logL = 0.;

  /* set up coeffs matrix */
  for ( i=0; i<Nmodels; i++ ){
    squared[i] = modelModel[i][i];
    coeffs[i][i] = -2.*dataModel[i];

    if( isinf(squared[i]) || isinf(coeffs[i][i]) ){ return -INFINITY; }

    for ( j=(i+1); j<Nmodels; j++ ){
      coeffs[i][j] = 2.*modelModel[i][j];
     
      if( isinf(coeffs[i][j]) ){ return -INFINITY; }
    }
  }

  for ( i=0; i<nm; i++ ){
    X = squared[i];
    invX = 1./X;
    invTwoX = 0.5*invX;
    invFourX = 0.25*invX;
    
    /* get the coefficients from the Y^2 term */
    for ( j=i; j<Nmodels; j++ ){
      for ( k=i; k<Nmodels; k++ ){
        /* add on new coefficients of squared terms */
        if ( j == i ){
          if ( k > j ){ squared[k] -= coeffs[j][k]*coeffs[j][k]*invFourX; }
          else if ( k == j ){ Z -= coeffs[j][k]*coeffs[j][k]*invFourX; }
        }
        else {
          if ( k == i ){ coeffs[j][j] -= coeffs[i][j]*coeffs[i][k]*invTwoX; }
          else if ( k > j ){ coeffs[j][k] -= coeffs[i][j]*coeffs[i][k]*invTwoX; }
        }
      }
    }
  }

  X = squared[nm];
  Y = coeffs[nm][nm];

  /* calculate analytic integral and get log likelihood */
  for ( i=0; i<Nmodels; i++ ){ logL -= 0.5*log(squared[i]); }

  logL -= 0.5*(Z - 0.25*Y*Y/X) / (sigma*sigma);

  logL += log(sigma);
  
  /* check whether final model is between 0 and infinity or -infinity and infinity */
  if (lastHalfRange == 1 ){
    logL += 0.5*(float)nm*AM_LN2PI + 0.5*AM_LNPI_2 + gsl_sf_log_erfc(0.5*Y/(sigma*sqrt(2.*X)));
  }
  else{ logL += 0.5*(float)Nmodels*AM_LN2PI; }

  return logL;
}


double log_marg_amp_full_C(int Nmodels, double modelModel[], double dataModel[], double sigma,
                           unsigned int lastHalfRange){
  /* coefficients of squares of model amplitudes */
  double squared[Nmodels];

  /* coefficients of model amplitudes */
  double coeffs[Nmodels][Nmodels];
  memset(coeffs, 0, sizeof(coeffs[0][0]) * Nmodels * Nmodels); /* initialise all values to zero */

  int i = 0, j = 0, k = 0, nm = Nmodels-1;

  double X = 0., invX = 0., invTwoX = 0., invFourX = 0., Y = 0., Z = 0., logL = 0.;

  /* set up coeffs matrix */
  for ( i=0; i<Nmodels; i++ ){
    squared[i] = modelModel[i*Nmodels + i];
    coeffs[i][i] = -2.*dataModel[i];

    if( isinf(squared[i]) || isinf(coeffs[i][i]) ){ return -INFINITY; }

    for ( j=(i+1); j<Nmodels; j++ ){
      coeffs[i][j] = 2.*modelModel[i*Nmodels + j];

      if( isinf(coeffs[i][j]) ){ return -INFINITY; }
    }
  }

  for ( i=0; i<nm; i++ ){
    X = squared[i];
    invX = 1./X;
    invTwoX = 0.5*invX;
    invFourX = 0.25*invX;

    /* get the coefficients from the Y^2 term */
    for ( j=i; j<Nmodels; j++ ){
      for ( k=i; k<Nmodels; k++ ){
        /* add on new coefficients of squared terms */
        if ( j == i ){
          if ( k > j ){ squared[k] -= coeffs[j][k]*coeffs[j][k]*invFourX; }
          else if ( k == j ){ Z -= coeffs[j][k]*coeffs[j][k]*invFourX; }
        }
        else {
          if ( k == i ){ coeffs[j][j] -= coeffs[i][j]*coeffs[i][k]*invTwoX; }
          else if ( k > j ){ coeffs[j][k] -= coeffs[i][j]*coeffs[i][k]*invTwoX; }
        }
      }
    }
  }

  X = squared[nm];
  Y = coeffs[nm][nm];

  /* calculate analytic integral and get log likelihood */
  for ( i=0; i<Nmodels; i++ ){ logL -= 0.5*log(squared[i]); }

  logL -= 0.5*(Z - 0.25*Y*Y/X) / (sigma*sigma);

  logL += log(sigma);

  /* check whether final model is between 0 and infinity or -infinity and infinity */
  if (lastHalfRange == 1 ){
    logL += 0.5*(double)nm*AM_LN2PI + 0.5*AM_LNPI_2 + gsl_sf_log_erfc(0.5*Y/(sigma*sqrt(2.*X)));
  }
  else{ logL += 0.5*(double)Nmodels*AM_LN2PI; }

  return logL;
}


/* For three amplitudes the marginalised likelihood ratio is given by:
 * \begin{equation}
 * \mathcal{L} = (2\pi)^{3/2}\sigma \frac{e^{-\frac{ d_h^2(f_g^2-FG) -2d_h(d_g f_g f_h - d_f f_h G -d_g F g_h + d_f f_g
g_h) + d_g^2(f_h^2 - FH) + 2d_f d_g(f_g H - f_h g_h) + d_f^2(g_h^2 - GH)}{2\sigma^2 (FGH + 2f_g f_h g_h -G f_h^2 
- H f_g^2 - F g_h^2) }}}{\sqrt{F}\sqrt{\frac{(FG - f_g^2)}{F}}\sqrt{\frac{FGH + 2f_g f_h g_h -G f_h^2 - F g_h^2 -
H f_g^2}{FG - f_g^2}}}
 * \end{equation}
 * for an integral between $[-\infty, \infty]$ for all amplitudes, and
   \begin{align}
   \mathcal{L} =& \sqrt{2}\pi^{3/2}\sigma \frac{e^{-\frac{ d_h^2(f_g^2-FG) -2d_h(d_g f_g f_h - d_f f_h G -d_g F g_h +
d_f f_g
g_h) + d_g^2(f_h^2 - FH) + 2d_f d_g(f_g H - f_h g_h) + d_f^2(g_h^2 - GH)}{2\sigma^2 (FGH + 2f_g f_h g_h -G f_h^2 
- H f_g^2 - F g_h^2) }}}{\sqrt{F}\sqrt{\frac{(FG - f_g^2)}{F}}\sqrt{\frac{FGH + 2f_g f_h g_h -G f_h^2 - F g_h^2 -
H f_g^2}{FG - f_g^2}}} \times \nonumber \\
 & \left( 1 + {\rm erf}\left( \frac{d_f f_g g_h + d_g f_g f_h + FG d_h - G d_f f_h - F d_g g_h - d_g
f_g^2}{\sqrt{2}\sigma(FG-f_g^2)\sqrt{\frac{FGH + 2f_g f_h g_h - F g_h^2 - G f_h^2 - h f_g^2}{FG-f_g^2}}} \right)
\right)
 \end{align}
 * for the final integral being between $[0, \infty]$. In these equations
 * we have used the following substitutions: $d_f = \sum_{i=1}^N d_i f_i$,
 * $d_g = \sum_{i=1}^N d_i g_i$, $d_h = \sum_{i=1}^N d_i h_i$, $f_g = \sum_{i=1}^N f_i g_i$
 * $f_h = \sum_{i=1}^N f_i h_i$, $g_h = \sum_{i=1}^N g_i h_i$, $F = \sum_{i=1}^N f_i^2$, 
 * $G = \sum_{i=1}^N g_i^2$, and $H = \sum_{i=1}^N H_i^2$ where $d$ is
 * the data and $f$, $g$ and $h$ are the three model components for which the
 * amplitudes have been marginalised over.  

/** This function explicitly writes out the marginalised likelihood ratio for a model consisting
 * of three components with the amplitudes marginalised over. This can be used to test the correctness
 * of the \c marginalise_amplitudes function */
double marginalise_three_amplitudes(double **modelModel, double *dataModel, double sigma, unsigned int lastHalfRange){
  double prefactor = 0.;
  double logL = 0., div = 0.;
  
  double F = modelModel[0][0], G = modelModel[1][1], H = modelModel[2][2];
  double fg = modelModel[0][1], fh = modelModel[0][2], gh = modelModel[1][2];
  double df = dataModel[0], dg = dataModel[1], dh = dataModel[2];
  
  prefactor = -0.5*log( F );
  prefactor -= 0.5*log( ( F*G - AM_SQUARE(fg) ) / F );
  prefactor -= 0.5*log( ( G*AM_SQUARE(fh) - 2.*fg*fh*gh + F*AM_SQUARE(gh) + H*AM_SQUARE(fg) - F*G*H ) / 
    ( AM_SQUARE(fg) - F*G ) );
  prefactor += log( sigma ) + 1.5*M_LNPI;
 
  if ( lastHalfRange == 1 ){
    double inerf = 0., X = ( AM_SQUARE(fg) - F*G );
    
    prefactor += 0.5*M_LN2;
    
    inerf = -dg*fg*fh + df*fh*G + dh*X + dg*F*gh - df*fg*gh;
    inerf /= ( sigma * M_SQRT2 * X * sqrt( ( AM_SQUARE(fh)*G - 2.*fg*fh*gh + 
                                             F*AM_SQUARE(gh) + H*AM_SQUARE(fg) - F*G*H ) / X ) ); 

    logL += log( 1. + gsl_sf_erf( inerf ) );
  }
  else{
    prefactor += 1.5*M_LN2;
  }
  
  div = 2.*AM_SQUARE(sigma) * ( G*AM_SQUARE(fh) - 2.*fg*fh*gh + H*AM_SQUARE(fg) + F*AM_SQUARE(gh) - F*G*H );
  
  logL += ( AM_SQUARE(dh)*(AM_SQUARE(fg)-F*G) - 2.*dh*(dg*fg*fh - df*fh*G - dg*F*gh + df*fg*gh) + 
            AM_SQUARE(dg)*(AM_SQUARE(fh)-F*H) + 2.*df*dg*(fg*H - fh*gh) + AM_SQUARE(df)*(AM_SQUARE(gh)-G*H) ) / div;

  logL += prefactor;
  
  return logL;
}


float marginalise_three_amplitudes_f(float **modelModel, float *dataModel, float sigma, unsigned int lastHalfRange){
  float prefactor = 0.;
  float logL = 0., div = 0.;
  
  float F = modelModel[0][0], G = modelModel[1][1], H = modelModel[2][2];
  float fg = modelModel[0][1], fh = modelModel[0][2], gh = modelModel[1][2];
  float df = dataModel[0], dg = dataModel[1], dh = dataModel[2];
  
  prefactor = -0.5*log( F );
  prefactor -= 0.5*log( ( F*G - AM_SQUARE(fg) ) / F );
  prefactor -= 0.5*log( ( G*AM_SQUARE(fh) - 2.*fg*fh*gh + F*AM_SQUARE(gh) + H*AM_SQUARE(fg) - F*G*H ) / 
    ( AM_SQUARE(fg) - F*G ) );
  prefactor += log( sigma ) + 1.5*M_LNPI;
 
  if ( lastHalfRange == 1 ){
    float inerf = 0., X = ( AM_SQUARE(fg) - F*G );
    
    prefactor += 0.5*M_LN2;
    
    inerf = -dg*fg*fh + df*fh*G + dh*X + dg*F*gh - df*fg*gh;
    inerf /= ( sigma * M_SQRT2 * X * sqrt( ( AM_SQUARE(fh)*G - 2.*fg*fh*gh + 
                                             F*AM_SQUARE(gh) + H*AM_SQUARE(fg) - F*G*H ) / X ) ); 

    logL += log( 1. + gsl_sf_erf( inerf ) );
  }
  else{
    prefactor += 1.5*M_LN2;
  }
  
  div = 2.*AM_SQUARE(sigma) * ( G*AM_SQUARE(fh) - 2.*fg*fh*gh + H*AM_SQUARE(fg) + F*AM_SQUARE(gh) - F*G*H );
  
  logL += ( AM_SQUARE(dh)*(AM_SQUARE(fg)-F*G) - 2.*dh*(dg*fg*fh - df*fh*G - dg*F*gh + df*fg*gh) + 
            AM_SQUARE(dg)*(AM_SQUARE(fh)-F*H) + 2.*df*dg*(fg*H - fh*gh) + AM_SQUARE(df)*(AM_SQUARE(gh)-G*H) ) / div;

  logL += prefactor;
  
  return logL;
}
