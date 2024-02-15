/* toeplitz.c */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <fftw3.h>

#ifdef __cplusplus
extern "C" {

int
fast_symm_toeplitz_prod (const size_t n,
			 const double * __restrict__ A1,
			 const double * __restrict__ x ,
			 double * __restrict__ b)
#else
int
fast_symm_toeplitz_prod (const size_t n,
			 const double * restrict A1,
			 const double * restrict x ,
			 double * restrict b)
#endif
{
  size_t k;   // padding
  size_t m;   /* the dimension of the augmented
                 matrix */
  
  /* Augmented arrays: A1 and x must be embedded
     into n+2*(k+1) arrays */
  double *mu;  // first row of augmented matrix
  double *y;   // augmented vector
  double *aux; // auxiliary array

  fftw_plan plan;

  /* The sizes needed to allocate memory for the
     augmented arrays */
  k = (n % 2 == 0) ? (n - 4) / 2 : (n - 3) / 2;
  m = n + 2 * (k + 1);

  /* Augmented arrays allocation */
  /*
  if ((mu  = (double *) malloc(m * sizeof(double))) == NULL)
      return 1;
  if ((y   = (double *) malloc(m * sizeof(double))) == NULL)
      return 2;

  if ((aux = (double *) malloc(m * sizeof(double))) == NULL)
      return 3;
  */
  
  if ((mu  = (double *) calloc(m, sizeof(double))) == NULL)
      return 1;
  if ((y   = (double *) calloc(m,  sizeof(double))) == NULL)
      return 2;

  if ((aux = (double *) calloc(m, sizeof(double))) == NULL)
      return 3;


  /* Construct of the the augmented arrays */
  /*
  memset((void *) (mu+n)   , 0, (m-n) * sizeof(double));
  memset((void *)  y       , 0, (k+1) * sizeof(double));
  memset((void *) (x+k+n+1), 0, (k+1) * sizeof(double));
  // TODO: INSTEAD THE ABOVE, TRY CALLOC
  // WARNING: THESE memsets are wrong !!!
  */
  
  memcpy((void *) mu, (const void *) A1, n * sizeof(double));
  for (size_t i = 2; i < n; ++i)
      mu[i-2] -= A1[i];
  memcpy((void *) (y+k+1), (const void *) x, n * sizeof(double));

  /* Compute the DST1 of mu and store in aux */
  plan = fftw_plan_r2r_1d(m, mu, aux,
    FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  /* Compute the DST1 of y and store in mu */
  plan = fftw_plan_r2r_1d(m, y, mu,
    FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  free(y); // y no more needed

  /* Multiply */
  {
    double scaling = 1.0 / ( 4.0 * (m + 1) );
    for (size_t i = 0; i < m; ++i)
        aux[i] *= (scaling * mu[i] / sin( (i+1) * M_PI / (m+1) ));
  }
  /* NOTE: need to multiply lamb by 1/2 rather
           sqrt(n+1)/2 because in fftw3 scaling
           of DST1 is 2 */

  /* Now compute the DST1 of xemb */
  plan = fftw_plan_r2r_1d(m, aux, mu,
    FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  free(aux); // aux no more needed

  /* Copy only necessary values to b */
  memcpy((void *) b, (const void *) (mu+k+1), n * sizeof(double));

  free(mu);

  return 0;
}

#ifdef __cplusplus
}  // extern "C" {
#endif
