/*   libb118
 *
 *   modules/linalg/include/b118/linalg/toeplitz/fstp_class.hpp
 *   
 *   Class for symmetric Toeplitz product
 *
 *   Copyright (C) 2024   Guilherme F. Fornel        <gffrnl@gmail.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <fftw3.h>
#include <cstddef>
#include <cstring>
#include <cmath>
#include <cassert>
#include <b118/utility.hpp>

struct fstp {
    double * aug_mat_first_row;
    double * aug_vec;
    double * aux;
    fftw_plan plan;
    std::size_t padding;
    std::size_t size;
    std::size_t aug_size;
    std::size_t capacity;

    void set_aug_size(std::size_t size) {
        padding  = (size % 2 == 0) ? (size - 4) / 2 : (size - 3) / 2;
        // (gffrnl) FFTW3 documentation says:
        //
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
        //     2^a 3^b 5^c 7^d 11^e 13^f, where e+f is either 0 or 1, and the
        //     other exponents are arbitrary. Other sizes are computed by means
        //     of a slow, general-purpose algorithm (which nevertheless retains
        //     O(n log n) performance even for prime sizes). (It is possible to
        //     customize FFTW for different array sizes; see Installation and
        //     Customization.)
        //     Transforms whose sizes are powers of 2 are especially fast.
        //
        //   - For a REDFT00 or RODFT00 transform kind in a dimension of size n,
        //     it is n-1 or n+1, respectively, that should be factorizable in
        //     the above form.
        //
        // ...
        //
        // (gffrnl) So, since we use RODFT00, we choose the augmented size
        // (gffrnl) of the form 2^a - 1:
        aug_size = b118::next_exp2(size + 2 * (padding + 1)) - 1;
    }

    void destroy_plan_and_free() {
        if (plan != nullptr) {
            fftw_destroy_plan(plan);
            plan = nullptr;
        }
        fftw_free(aux);
        aux = nullptr;
        fftw_free(aug_vec);
        aug_vec = nullptr;
        fftw_free(aug_mat_first_row);
        aug_mat_first_row = nullptr;
    }

    void alloc_fftw(std::size_t n) {
        aug_mat_first_row =
            static_cast<double *>(fftw_malloc(n * sizeof(double)));
        aug_vec =
            static_cast<double *>(fftw_malloc(n * sizeof(double)));
        aux =
            static_cast<double *>(fftw_malloc(n * sizeof(double)));
    }

 public:
    explicit fstp(std::size_t n) : aug_mat_first_row(nullptr),
                                   aug_vec(nullptr),
                                   aux(nullptr),
                                   plan(nullptr),
                                   padding(0),
                                   size(n),
                                   aug_size(0),
                                   capacity(n) {
        set_aug_size(size);
        if (capacity < aug_size) {
            capacity = aug_size;
        }
        alloc_fftw(capacity);
        // size the plan
        plan = fftw_plan_r2r_1d(aug_size, aug_mat_first_row, aug_vec,
                                FFTW_RODFT00, FFTW_ESTIMATE);
    }

    ~fstp() {
        destroy_plan_and_free();
    }

    std::size_t get_size() const { return size; }

    fstp& resize(std::size_t size) {
        if (this->size != size) {
            set_aug_size(size);
            if (aug_size > capacity) {
                destroy_plan_and_free();
                alloc_fftw(aug_size);
                capacity = aug_size;
            }
            // size the plan
            plan = fftw_plan_r2r_1d(aug_size, aug_mat_first_row, aug_vec,
                                    FFTW_RODFT00, FFTW_ESTIMATE);
            this->size = size;
        }
        return *this;
    }

    fstp& multiply(double const * const __restrict__ mat_first_row,
                   double const * const __restrict__ vec,
                   double       * const __restrict__ prod) {
        assert(plan != nullptr);

        // Construct of the the augmented arrays
        std::memcpy(static_cast<void *>(aug_mat_first_row),
                    static_cast<void const *>(mat_first_row),
                    size * sizeof(double));
        for (std::size_t i = 2; i < size; ++i)
            aug_mat_first_row[i-2] -= mat_first_row[i];
        std::memcpy(static_cast<void *>(aug_vec + padding + 1),
                    static_cast<void const *>(vec),
                    size * sizeof(double));

        // Populate with 0 the rest of
        std::memset(static_cast<void *>(aug_mat_first_row + size),
                    0,
                    (aug_size - size) * sizeof(double));
        std::memset(static_cast<void *>(aug_vec),
                    0,
                    (padding + 1) * sizeof(double));
        std::memset(static_cast<void *>(aug_vec + padding + 1 + size),
                    0,
                    (aug_size - (padding + 1 + size)) * sizeof(double));

        // Compute the DST1 of aug_mat_first_row and store in aux
        fftw_execute_r2r(plan, aug_mat_first_row, aux);

        // Compute the DST1 of y and store in aug_mat_first_row
        fftw_execute_r2r(plan, aug_vec, aug_mat_first_row);

        // Multiply
        {
            double const scaling = 1.0 / (4.0 * (aug_size + 1));
            double const delta_theta = M_PI / (aug_size + 1);
            for (std::size_t i = 0; i < aug_size; ++i)
                aux[i] *= (scaling * aug_mat_first_row[i]
                        / std::sin((i + 1) * delta_theta));
        }
        // NOTE: TODO(gffrnl) REWRITE THIS NOTE
        // need to multiply lamb by 1/2 rather
        //          sqrt(n+1)/2 because in fftw3 scaling
        //          of DST1 is 2

        // Now compute the DST1 of aux
        fftw_execute_r2r(plan, aux, aug_mat_first_row);

        // Copy only necessary values to prod
        std::memcpy(static_cast<void *>(prod),
                    static_cast<void const *>(aug_mat_first_row + padding + 1),
                    size * sizeof(double));

        return *this;
    }

    void cleanup() { fftw_cleanup(); }
};
