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

#include <ampmarginaliser.h>

#define TESTLEN 50
#define TESTMODEL1(x) x
#define TESTMODEL2(x) x*x
#define TESTMODEL3(x) exp(-(x*x)/2.)
#define TESTMODEL4(x) sin(x)
#define TESTFRACERR 1e-9

int test_marginalisation();
int test_marginalisation_except_final();

/** Test the marginalisation performed by \c marginalise_amplitudes against that of \c marginalise_three_amplitudes
 */
int test_marginalisation(){
  double data[TESTLEN];
  double sigma = 2.0; /* standard deviation of the data */
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


int test_marginalisation_except_final(){
  double data[TESTLEN];
  double sigma = 2.0; /* standard deviation of the data */
  double model[4][TESTLEN];
  double amplitudes[4] = {1., 2., 3., 4.};
  double post1 = 0., post2 = 0.; /* posterior samples */
  double amp1 = 0., amp2 = 0.;

  const gsl_rng_type * T;
  gsl_rng *rng;

  int i = 0, j = 0, k = 0, n = 0;

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

  double *dataModel = calloc(4, sizeof(double));
  double **modelModel = NULL;

  /* create data and models (just random Gaussian noise) */
  for ( i=0; i<TESTLEN; i++ ){
    data[i] = gsl_ran_gaussian(rng, sigma);
    model[0][i] = TESTMODEL1((double)i/(double)(TESTLEN-1.));
    model[1][i] = TESTMODEL2((double)i/(double)(TESTLEN-1.));
    model[2][i] = TESTMODEL3((double)i/(double)(TESTLEN-1.));
    model[3][i] = TESTMODEL4((double)i/(double)(TESTLEN-1.));
  }

  gsl_rng_free(rng);

  modelModel = calloc(4, sizeof(double*));

  /* create dataModel and modelModel arrays */
  for ( i=0; i<4; i++ ){
    modelModel[i] = calloc(4, sizeof(double));
  }

  for ( n=0; n<4; n++ ){
    for ( i=0; i<4; i++ ){
      /* create dataModel cross terms */
      if ( i == 3 ){ amp1 = amplitudes[n]; }
      else { amp1 = 1.; }

      for ( k=0; k<TESTLEN; k++ ){ dataModel[i] += data[k]*amp1*model[i][k]; }

      for ( j=i; j<4; j++ ){
        if ( j == 3 ){ amp2 = amplitudes[n]; }
        else { amp2 = 1.; }

        for ( k=0; k<TESTLEN; k++ ){ modelModel[i][j] += amp1*model[i][k]*amp2*model[j][k]; }
      }
    }

    /* get log likelihood ratio from marginalise_amplitudes_except_final */
    post1 = marginalise_amplitudes_except_final(4, modelModel, dataModel, sigma);

    /* get log likelihood ratio from marginalise_three_amplitudes */
    post2 = marginalise_three_amplitudes_exclude_final(modelModel, dataModel, sigma);

    fprintf(stdout, "post1[%d] (all integrated from -inf to inf) marginalise_amplitudes_except_final():\t%.9lf\n", n, post1);
    fprintf(stdout, "post2[%d] (all integrated from -inf to inf) marginalise_three_amplitudes_except_final():\t%.9lf\n", n, post2);

    /* test whether the values are the same to within a tolerable fractional error */
    if ( fabs((post1-post2)/post1) > TESTFRACERR ){
      fprintf(stderr, "Probabilities are not the same!\n");
      free(dataModel);
      for (i=0; i<4; i++){ free(modelModel[i]); }
      free(modelModel);
      return 0;
    }
  }

  free(dataModel);
  for (i=0; i<4; i++){ free(modelModel[i]); }
  free(modelModel);

  return 1; /* successful comparison */
}


int main(void){
  int i = test_marginalisation();

  if ( i ){
    fprintf(stdout, "SUCCESS... the likelihood ratios are the same! :D\n");
  }
  else { return -1; }

  i = test_marginalisation_except_final();

  if ( i ){
    fprintf(stdout, "SUCCESS... the probabilities are the same! :D\n");
  }
  else { return -1; }

  return 0;
}

