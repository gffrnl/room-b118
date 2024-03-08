/*   libb118
 *
 *   modules/linalg/include/b118/linalg/toeplitz/fast_symm_toeplitz_product.hpp
 *   
 *   Fast symmetric Toeplitz matrix-vector product
 *
 *   Copyright (C) 2024   Guilherme F. Fornel        <gffrnl@gmail.com>
 *                        Fabio Souto de Azevedo     <fazedo@gmail.com>
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

#include <cmath>
#include <algorithm>
#include <b118/numbers.hpp>
#include "./fsymmtp_manager.hpp"

template<typename Real>
struct fftw_r2r_wrapper final{
    static void execute(b118::fftw::plan<Real> const & plan,
                        Real * in, Real * out) {
        fftwq_execute_r2r(plan.raw, in, out);
    }
};
template<>
struct fftw_r2r_wrapper<float> final {
    static void execute(b118::fftw::plan<float> const & plan,
                        float * in, float * out) {
        fftwf_execute_r2r(plan.raw, in, out);
    }
};
template<>
struct fftw_r2r_wrapper<double> final {
    static void execute(b118::fftw::plan<double> const & plan,
                        double * in, double * out) {
        fftw_execute_r2r(plan.raw, in, out);
    }
};
template<>
struct fftw_r2r_wrapper<long double> final {
    static void execute(b118::fftw::plan<long double> const & plan,
                        long double * in, long double * out) {
        fftwl_execute_r2r(plan.raw, in, out);
    }
};

template<
    typename Real,
    class InputBidirIt1, class InputBidirIt2,
    class OutputBidirIt
>
OutputBidirIt fast_symm_toeplitz_product(
    InputBidirIt1 row_beg, InputBidirIt1 row_end,
    InputBidirIt2 vec_beg,
    OutputBidirIt prod_beg, bool inplace = false) {
    auto const len = std::distance(row_beg, row_end);

    {
        if (len < static_cast<decltype(len)>(0))
            throw std::invalid_argument("fast_symm_toeplitz_product(): "
                                        "std::distance(row_beg, row_end) < 0");

        if (len == static_cast<decltype(len)>(0))
            return prod_beg;

        if (len == static_cast<decltype(len)>(1)) {
            *prod_beg = (*row_beg) * (*vec_beg);
            return std::next(prod_beg);
        }

        fsymmtp_manager<Real>::resize(static_cast<std::size_t>(len));
    }

    std::size_t size     = fsymmtp_manager<Real>::size();
    std::size_t padding  = fsymmtp_manager<Real>::padding();
    std::size_t aug_size = fsymmtp_manager<Real>::aug_size();

    Real * const aug_row = fsymmtp_manager<Real>::aug_row();
    Real * const aug_vec = fsymmtp_manager<Real>::aug_vec();
    Real * const aux     = fsymmtp_manager<Real>::aux();

    b118::fftw::plan<Real> const plan = fsymmtp_manager<Real>::plan();

    // Construct the augmented arrays
    std::copy(row_beg, row_end, aug_row);
    for (std::size_t i = 2; i < size; ++i)
        aug_row[i-2] -= (*std::next(row_beg, i));
    std::copy_n(vec_beg, size, aug_vec + padding + 1);

    // Populate with 0 the rest of
    std::fill_n(aug_row + size, aug_size - size, 0);
    std::fill_n(aug_vec, padding + 1, 0);
    std::fill_n(aug_vec + padding + 1 + size,
                aug_size - (padding + 1 + size), 0);


    // // Populate with 0 the rest of
    // std::memset(static_cast<void *>(aug_mat_first_row + size),
    //             0,
    //             (aug_size - size) * sizeof(double));
    // std::memset(static_cast<void *>(aug_vec),
    //             0,
    //             (padding + 1) * sizeof(double));
    // std::memset(static_cast<void *>(aug_vec + padding + 1 + size),
    //             0,
    //             (aug_size - (padding + 1 + size)) * sizeof(double));

    // Compute the DST1 of aug_mat_first_row and store in aux
    // fftw_execute_r2r(plan.ptr, aug_row, aux);
    fftw_r2r_wrapper<Real>::execute(plan, aug_row, aux);

    // Compute the DST1 of y and store in aug_mat_first_row
    // fftw_execute_r2r(plan.ptr, aug_vec, aug_row);
    fftw_r2r_wrapper<Real>::execute(plan, aug_vec, aug_row);

    // Multiply
    {
        Real const scaling = static_cast<Real>(1)
            / static_cast<Real>(4 * (aug_size + 1));
        Real const delta_theta = b118::numbers::pi_v<Real>
            / static_cast<Real>(aug_size + 1);
        for (std::size_t i = 0; i < aug_size; ++i)
            aux[i] *= (scaling * aug_row[i]
                    / sin(static_cast<Real>((i + 1) * delta_theta)));
    }
    // NOTE: TODO(gffrnl) REWRITE THIS NOTE
    // need to multiply lamb by 1/2 rather
    //          sqrt(n+1)/2 because in fftw3 scaling
    //          of DST1 is 2

    // Now compute the DST1 of aux
    // fftw_execute_r2r(plan.ptr, aux, aug_row);
    fftw_r2r_wrapper<Real>::execute(plan, aux, aug_row);

    // Copy only necessary values to prod
    if (!inplace) {
        std::copy_n(aug_row + padding + 1, size, prod_beg);
    } else {
        for (std::size_t i = 0; i < size; ++i) {
            // b[i] += mu[i + k + 1];
            *(prod_beg + i) += *(aug_row + i +padding + 1);
        }
    }

    return std::next(prod_beg, size);
}
