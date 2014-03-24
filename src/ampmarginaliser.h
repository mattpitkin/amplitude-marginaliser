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

#include <math.h>
#include <string.h>

/* GSL headers */
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define AM_SQUARE(x) (x*x)
#define AM_LN2PI (M_LN2 + M_LNPI)
#define AM_LNPI_2 (M_LNPI - M_LN2)

/* define functions */
double marginalise_amplitudes(int Nmodels, double **modelModel, double *dataModel, double sigma,
                              unsigned int lastHalfRange);
float marginalise_amplitudes_f(int Nmodels, float **modelModel, float *dataModel, float sigma,
                               unsigned int lastHalfRange);

/* version of marginalise amplitudes where the modelModel array is held in a 1D vector rather than a 2D array */
double marginalise_amplitudes_linear(int Nmodels, double modelModel[], double dataModel[], double sigma,
                                     unsigned int lastHalfRange);

double marginalise_three_amplitudes(double **modelModel, double *dataModel, double sigma, unsigned int lastHalfRange);
float marginalise_three_amplitudes_f(float **modelModel, float *dataModel, float sigma, unsigned int lastHalfRange);
