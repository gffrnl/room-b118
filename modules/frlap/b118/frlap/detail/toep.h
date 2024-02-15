/* toeplitz.h */
#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void symm_toeplitz_prod(size_t const                      n,
                        double const * const __restrict__ A1,
                        double const * const __restrict__ x,
                        double       * const __restrict__ b);

int fast_symm_toeplitz_prod(size_t const                      n,
                            double const * const __restrict__ A1,
                            double const * const __restrict__ x,
                            double       * const __restrict__ b);


void durbin(size_t const                      n,
            double const * const __restrict__ r,
            double       * const __restrict__ y);

int levinson(size_t const                      n,
             double const * const __restrict__ A1,
             double       * const __restrict__ x,
             double const * const __restrict__ b);

int levinson_in_place(size_t const                      n,
                      double const * const __restrict__ A1,
                      double       * const __restrict__ x,
                      double const * const __restrict__ b);

#ifdef __cplusplus
}  // end extern "C"
#endif
