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

/** \file ampmarginaliser.c
 */

#include "ampmarginaliser.h"

double marginalise_amplitudes(int Nmodels, double **modelModel, double *dataModel, double sigma,
                              unsigned int lastHalfRange){
  /**  \var Coefficients of squares of model amplitudes */
  double squared[Nmodels];

  /** \var coeff An array to hold the coefficients of the model amplitudes */
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
  
  /** For e.g. five models components the \c coeff matrix would be equivalent to:
   \f[ \left(
   \begin{array}{ccccc}
     -2DF & 2fg   & 2fh  &  2fi &  2fj \\
      0     & -2DG  & 2gh  &  2gi &  2gj \\
      0     & 0     & -2DH &  2hi &  2hj \\
      0     & 0     & 0    & -2DI &  2ij \\
      0     & 0     & 0    & 0    & -2DJ
   \end{array}
   \right)
   \f]
   * where e.g. \f$DF = \sum d_if_i\f$ and \f$fg = \sum f_ig_i\f$, and where for
   * models with amplitudes \f$A\f$, \f$B\f$, \f$C\f$, \f$D\f$ and \f$E\f$ these
   * are the coefficients from the likelihood for each, with rows corresponding to
   * coefficients of each amplitude \f$A\f$ through \f$E\f$, and columns corresponding
   * to cross term amplitudes, e.g.:
   * <br><tt>coeffs[0][0]</tt> is terms with coefficients just of \f$A\f$,
   * <br><tt>coeffs[0][1]</tt> is terms with coefficients of \f$A\f$ and \f$B\f$,
   * <br><tt>coeffs[0][1]</tt> is terms with coefficients of \f$A\f$ and \f$C\f$
   * <br><tt>...</tt>
   * <br><tt>coeffs[2][2]</tt> is terms with coefficients of just \f$C\f$,
   * <br><tt>coeffs[2][3]</tt> is terms with coefficients of \f$C\f$ and \f$D\f$,
   * <br><tt>...</tt>
   * <br><tt>coeffs[4][4]</tt> is terms with coefficients just of \f$E\f$.
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


double marginalise_amplitudes_linear(int Nmodels, double modelModel[], double dataModel[], double sigma,
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
