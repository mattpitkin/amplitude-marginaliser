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

#include <ampmarginaliser/ampmarginaliser.h>

#define TESTLEN 50
#define TESTMODEL1(x) x
#define TESTMODEL2(x) x*x
#define TESTMODEL3(x) exp(-(x*x)/2.)
#define TESTFRACERR 1e-9

int test_marginalisation();

/** Test the marginalisation performed by \c marginalise_amplitudes against that of \c marginalise_three_amplitudes
 */
int test_marginalisation(){
  double data[TESTLEN];
  double sigma = 1.0; /* standard deviation of the data */
  double model[3][TESTLEN];
  
  const gsl_rng_type * T;
  gsl_rng *rng;
 
  int i = 0, j = 0, k = 0;

  double logL1, logL2;
  
  /* get random seed */
  FILE *devrandom = NULL;
  struct timeval tv;
  int randomseed;
  if ( (devrandom = fopen("/dev/random","r")) == NULL ) {
    gettimeofday( &tv, 0 );
    randomseed = tv.tv_sec + tv.tv_usec;
  }
  else {
    if( fread(&randomseed, sizeof(randomseed), 1, devrandom) != 1 ){
      fprintf(stderr, "Error... could not read random seed\n");
      return 0;
    }
    fclose( devrandom );
  }
  
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  gsl_rng_set(rng, randomseed);
  
  double *dataModel = calloc(3, sizeof(double));
  double **modelModel = NULL;

  /* create data and models (just random Gaussian noise) */
  for ( i=0; i<TESTLEN; i++ ){
    data[i] = gsl_ran_gaussian(rng, sigma);
    model[0][i] = TESTMODEL1((double)i/(double)(TESTLEN-1.));
    model[1][i] = TESTMODEL2((double)i/(double)(TESTLEN-1.));
    model[2][i] = TESTMODEL3((double)i/(double)(TESTLEN-1.));
  }
  
  gsl_rng_free(rng);
  
  modelModel = calloc(3, sizeof(double*));
  
  /* create dataModel and modelModel arrays */
  for ( i=0; i<3; i++ ){
    modelModel[i] = calloc(3, sizeof(double));
    
    /* create dataModel cross terms */
    for ( k=0; k<TESTLEN; k++ ){ dataModel[i] += data[k]*model[i][k]; }
    
    for ( j=i; j<3; j++ ){
      for ( k=0; k<TESTLEN; k++ ){ modelModel[i][j] += model[i][k]*model[j][k]; }
    }
  }
  
  /* get log likelihood ratio from marginalise_amplitudes, with lastHalfRange = 0 */
  logL1 = marginalise_amplitudes(3, modelModel, dataModel, sigma, 0);
  
  /* get log likelihood ratio from marginalise_three_amplitudes, with lastHalfRange = 0 */
  logL2 = marginalise_three_amplitudes(modelModel, dataModel, sigma, 0);
  
  fprintf(stdout, "logL (all integrated from -inf to inf) marginalise_amplitudes():\t%.9lf\n", logL1);
  fprintf(stdout, "logL (all integrated from -inf to inf) marginalise_three_amplitudes():\t%.9lf\n", logL2);
  
  /* test whether the values are the same to within a tolerable fractional error */
  if ( fabs((logL1-logL2)/logL1) > TESTFRACERR ){
    fprintf(stderr, "Likelihoods are not the same!\n");
    free(dataModel);
    for (i=0; i<3; i++){ free(modelModel[i]); }
    free(modelModel);
    return 0;
  }
  
  /* get log likelihood ratio from marginalise_amplitudes, with lastHalfRange = 0 */
  logL1 = marginalise_amplitudes(3, modelModel, dataModel, sigma, 1);
  
  /* get log likelihood ratio from marginalise_three_amplitudes, with lastHalfRange = 0 */
  logL2 = marginalise_three_amplitudes(modelModel, dataModel, sigma, 1);
  
  fprintf(stdout, "logL (last integrated from 0 to inf) marginalise_amplitudes():\t\t%.9lf\n", logL1);
  fprintf(stdout, "logL (last integrated from 0 to inf) marginalise_three_amplitudes():\t%.9lf\n", logL2);

  /* test whether the values are the same to within a tolerable fractional error */
  if ( fabs((logL1-logL2)/logL1) > TESTFRACERR ){
    fprintf(stderr, "Likelihoods are not the same!\n");
    free(dataModel);
    for (i=0; i<3; i++){ free(modelModel[i]); }
    free(modelModel);
    return 0;
  }
  
  free(dataModel);
  for (i=0; i<3; i++){ free(modelModel[i]); }
  free(modelModel);
  
  return 1; /* successful comparison */
}


int main(void){
  int i = test_marginalisation();
  
  if ( i ){ 
    fprintf(stdout, "SUCCESS... the likelihood ratios are the same! :D\n");
    return 0; }
  else { return -1; }
}
