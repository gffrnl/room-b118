// Copyright 2024 Guilherme F. Fornel

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <float.h>


#include <fftw3.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cmath>

void fast_symm_toeplitz_prod(std::size_t const                      n,
                             double      const * const __restrict__ A1,
                             double      const * const __restrict__ x,
                             double            * const __restrict__ b) {
  std::size_t k;  // padding
  std::size_t m;  // the dimension of the augmented matrix

  // Augmented arrays: A1 and x must be embedded into n+2*(k+1) arrays
  double * mu;   // first row of augmented matrix
  double * y;    // augmented vector
  double * aux;  // auxiliary array

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
      abort();
  if ((y   = (double *) calloc(m,  sizeof(double))) == NULL)
      abort();

  if ((aux = (double *) calloc(m, sizeof(double))) == NULL)
      abort();


  /* Construct of the the augmented arrays */
  /*
  memset((void *) (mu+n)   , 0, (m-n) * sizeof(double));
  memset((void *)  y       , 0, (k+1) * sizeof(double));
  memset((void *) (x+k+n+1), 0, (k+1) * sizeof(double));
  // TODO: INSTEAD THE ABOVE, TRY CALLOC
  // WARNING: THESE memsets are wrong !!!
  */

  std::memcpy(static_cast<void *>(mu),    static_cast<void const *>(A1),
              n * sizeof(double));
  for (size_t i = 2; i < n; ++i)
      mu[i-2] -= A1[i];
  std::memcpy(static_cast<void *>(y+k+1), static_cast<void const *>(x),
              n * sizeof(double));

  // Compute the DST1 of mu and store in aux
  plan = fftw_plan_r2r_1d(m, mu, aux, FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  // fftw_cleanup();

  // Compute the DST1 of y and store in mu
  plan = fftw_plan_r2r_1d(m, y, mu, FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  // fftw_cleanup();

  std::free(y);  // y no more needed

  // Multiply
  {
    double const scaling = 1.0 / (4.0 * (m + 1));
    double const delta_theta = M_PI / (m + 1);
    for (std::size_t i = 0; i < m; ++i)
        aux[i] *= (scaling * mu[i] / sin((i + 1) * delta_theta));
  }
  // NOTE: TODO(grffrnl) REWRITE THIS NOTE
  // need to multiply lamb by 1/2 rather
  //          sqrt(n+1)/2 because in fftw3 scaling
  //          of DST1 is 2

  /* Now compute the DST1 of xemb */
  plan = fftw_plan_r2r_1d(m, aux, mu, FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  std::free(aux);  // aux no more needed

  // Copy only necessary values to b
  std::memcpy(static_cast<void *>(b), static_cast<void const *>(mu+k+1),
              n * sizeof(double));

  std::free(mu);
}
