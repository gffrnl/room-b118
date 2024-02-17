// Copyright 2024 Guilherme F. Fornel

#include <stdio.h>
#include <stdint.h>
#include <float.h>


#include <fftw3.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <b118/utility.hpp>

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

    // The sizes needed to allocate memory for the
    //    augmented arrays
    k = (n % 2 == 0) ? (n - 4) / 2 : (n - 3) / 2;
    m = n + 2 * (k + 1);

    // https://www.fftw.org/fftw3_doc/
    //
    //
    // 4.3.5 Real-to-Real Transforms
    //
    // ...
    //
    // fftw_plan fftw_plan_r2r_1d(int n, double *in, double *out,
    //                            fftw_r2r_kind kind, unsigned flags);
    //
    // * n gives the (physical) size of the transform dimensions.
    //   They can be any positive integer.
    //
    //   - FFTW is generally best at handling sizes of the form
    //     2a 3b 5c 7d 11e 13f, where e+f is either 0 or 1, and the other
    //     exponents are arbitrary. Other sizes are computed by means of a slow,
    //     general-purpose algorithm (which nevertheless retains O(n log n)
    //     performance even for prime sizes). (It is possible to customize FFTW
    //     for different array sizes; see Installation and Customization.)
    //     Transforms whose sizes are powers of 2 are especially fast.
    //
    //   - For a REDFT00 or RODFT00 transform kind in a dimension of size n,
    //     it is n-1 or n+1, respectively, that should be factorizable in
    //     the above form.
    //
    // ...
    //
    std::size_t const m2 = b118::next_exp2(m) - 1;

    mu  = static_cast<double *>(fftw_malloc(sizeof(double) * m2));
    y   = static_cast<double *>(fftw_malloc(sizeof(double) * m2));
    aux = static_cast<double *>(fftw_malloc(sizeof(double) * m2));

    // Construct of the the augmented arrays
    std::memcpy(static_cast<void *>(mu),    static_cast<void const *>(A1),
                n * sizeof(double));
    for (std::size_t i = 2; i < n; ++i)
        mu[i-2] -= A1[i];
    std::memcpy(static_cast<void *>(y + k + 1), static_cast<void const *>(x),
                n * sizeof(double));

    // Populate with 0 the rest of
    std::memset(static_cast<void *>(mu + n),
                0,
                (m2 - n) * sizeof(double));
    std::memset(static_cast<void *>(y),
                0,
                (k + 1) * sizeof(double));
    std::memset(static_cast<void *>(y + k + 1 + n),
                0,
                (m2 - (k + 1 + n)) * sizeof(double));

    // Compute the DST1 of mu and store in aux
    plan = fftw_plan_r2r_1d(m2, mu, aux, FFTW_RODFT00, FFTW_ESTIMATE);
    fftw_execute_r2r(plan, mu, aux);

    // Compute the DST1 of y and store in mu
    fftw_execute_r2r(plan, y, mu);

    fftw_free(y);

    // Multiply
    {
        double const scaling = 1.0 / (4.0 * (m2 + 1));
        double const delta_theta = M_PI / (m2 + 1);
        for (std::size_t i = 0; i < m2; ++i)
            aux[i] *= (scaling * mu[i] / sin((i + 1) * delta_theta));
    }
    // NOTE: TODO(grffrnl) REWRITE THIS NOTE
    // need to multiply lamb by 1/2 rather
    //          sqrt(n+1)/2 because in fftw3 scaling
    //          of DST1 is 2

    // Now compute the DST1 of aux
    fftw_execute_r2r(plan, aux, mu);
    fftw_destroy_plan(plan);

    fftw_free(aux);

    // Copy only necessary values to b
    std::memcpy(static_cast<void *>(b), static_cast<void const *>(mu + k + 1),
                n * sizeof(double));

    fftw_free(mu);

    fftw_cleanup();
}
