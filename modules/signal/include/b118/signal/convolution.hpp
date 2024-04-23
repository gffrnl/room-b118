/*   libb118
 *
 *   modules/signal/include/b118/signal/convolution.hpp
 *   
 *   COnvolution and Cross-correlation
 *
 *   Copyright (C) 2024   Fabio Souto de Azevedo     <fazedo@gmail.com>
 *                        Guilherme F. Fornel        <gffrnl@gmail.com>
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

// Copyright 2024 Fabio Souto de Azevedo <fazedo@gmail.com>

#pragma once

// #include <x86_64-linux-gnu/cblas.h>
#include <cblas.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <typeinfo>
#include <cassert>
#include <algorithm>
#include <execution>
#include <b118/fftw.hpp>

namespace b118 { namespace signal {

// The input is y and mu, and the output is frLap_y.
// mu must have at least na + n0 - 1 elements
// y must have at least na elements
// The convolution is given by:
// frLap_y[i] = sum(mu[(ja + i - j] * y[j], j = 0 ... na - 1), i = 0 .. n0 - 1
// The cross-correlation is given by:
// frLap_y[i] = sum(mu[n0 - 1 - i + j] * y[j + jb], j = 0 ... nb - 1),
//                                                  i = 0 ... n0 - 1

template<typename Real>
class convolution final : private b118::fftw::managed<Real> {
    Real * in;
    Real * kernel;
    Real * out;
    b118::fftw::complex<Real> * in_fft;
    b118::fftw::complex<Real> * kernel_fft;
    b118::fftw::plan<Real, b118::fftw::r2c> plan_r2c;
    b118::fftw::plan<Real, b118::fftw::c2r> plan_c2r;
    std::size_t conv_size = 0;  // The number of real samples in current plan
    std::size_t capacity = 0;

    void elementwise_product(bool conj = false) {
        if (conj == false) {
            for (std::size_t i = 0; i < conv_size/2 + 1; ++i) {
                Real re = + in_fft[i][0] * kernel_fft[i][0]
                          - in_fft[i][1] * kernel_fft[i][1];
                Real im = + in_fft[i][0] * kernel_fft[i][1]
                          + in_fft[i][1] * kernel_fft[i][0];
                in_fft[i][0] = re;
                in_fft[i][1] = im;
            }
        } else {
            for (std::size_t i = 0; i < conv_size/2 + 1; ++i) {
                Real re = + in_fft[i][0] * kernel_fft[i][0]
                          + in_fft[i][1] * kernel_fft[i][1];
                Real im = - in_fft[i][0] * kernel_fft[i][1]
                          + in_fft[i][1] * kernel_fft[i][0];
                in_fft[i][0] = re;
                in_fft[i][1] = im;
            }
        }
    }

 public:
    explicit convolution(std::size_t N) : capacity(N) {
        // N is the logical size, i.e. the number of real samples
        // N must be at least na + n0 - 1
        // We may initialize with a large N and then
        // redo the plans with the actual size of the convolution.
        // This avoids the overhead of allocation and reallocating memory.
        in         = b118::fftw::alloc<Real>::real(N);
        kernel     = b118::fftw::alloc<Real>::real(N);
        out        = b118::fftw::alloc<Real>::real(N);
        in_fft     = b118::fftw::alloc<Real>::complex(N/2 + 1);
        kernel_fft = b118::fftw::alloc<Real>::complex(N/2 + 1);
    }

    ~convolution() {
        // b118::fftw::destroy_plan(plan_r2c);
        // b118::fftw::destroy_plan(plan_c2r);
        b118::fftw::free<Real>::real(in);
        b118::fftw::free<Real>::real(kernel);
        b118::fftw::free<Real>::real(out);
        b118::fftw::free<Real>::complex(in_fft);
        b118::fftw::free<Real>::complex(kernel_fft);
    }

    // convolution& create_plans(std::size_t conv_size) {
    //     if (plan_r2c.raw != nullptr) {
    //         b118::fftw::destroy_plan(plan_r2c);
    //     }
    //     if (plan_c2r.raw != nullptr) {
    //         b118::fftw::destroy_plan(plan_c2r);
    //     }
    //     // conv_size the number of real samples.
    //     plan_r2c = b118::fftw::dft::r2c::create_plan(
    //         conv_size, in, in_fft, FFTW_ESTIMATE);
    //     plan_c2r = b118::fftw::dft::c2r::create_plan(
    //         conv_size, in_fft, out, FFTW_ESTIMATE);
    //     this->conv_size = conv_size;
    //     return *this;
    // }

    convolution& create_plans(std::size_t conv_size) {
        // conv_size the number of real samples.
        plan_r2c = b118::fftw::plan<Real, b118::fftw::r2c>(
            conv_size, in, in_fft, FFTW_ESTIMATE);        
        plan_c2r = b118::fftw::plan<Real, b118::fftw::c2r>(
            conv_size, in_fft, out, FFTW_ESTIMATE);
        this->conv_size = conv_size;
        return *this;
    }

    convolution& conv(std::size_t output_size,
                      std::size_t input_size,
                      Real const * mu,
                      Real const * input,
                      Real *output,
                      bool is_cross_correlation = false) {

        // Copy y to in
        if (is_cross_correlation == true) {
            std::fill(in, in + conv_size - input_size, 0.0);
            std::copy(input, input + input_size, in + conv_size - input_size);
        } else {
            std::copy(input, input + input_size, in);
            std::fill(in + input_size, in + conv_size, 0.0);
        }

        // mu to kernel. Note that we are copying na + n0 - 1 elements
        // Hence the kernel is always larger than y, which is not a problem
        // at all, but it may be a bit confusing.
        // Warning: we need to copy only the slice of mu that we are going to use
        // The indices must be calculated before calling this function
        auto ctn = std::copy(mu,
                             mu + input_size + output_size - 1,
                             kernel);

        std::fill(ctn, kernel + conv_size, 0.0);

        // assert(conv_size <= max_size);
        // Execute plan_r2c for each
        // b118::fftw::dft::r2c::execute(plan_r2c, in, in_fft);
        // b118::fftw::dft::r2c::execute(plan_r2c, kernel, kernel_fft);
        plan_r2c.execute(in, in_fft);
        plan_r2c.execute(kernel, kernel_fft);

        // Elementwise product
        elementwise_product(is_cross_correlation);

        // Execute plan_c2r, i.e., inverse transform
        // From doc: the c2r transform destroys its input array even for out-of-place transforms.
        // b118::fftw::dft::c2r::execute(plan_c2r, in_fft, out);
        plan_c2r.execute(in_fft, out);

        // Copy and scale output from correctio positions according to the operation
        Real* beginning_valid_window =
                        out + (is_cross_correlation ? 0 : input_size - 1);

        // Blas DAXPY function computes a constant times a vector plus a vector.
        // output <- output + beginning_valid_window/conv_size

        std::transform(std::execution::par,
                       const_cast<Real *>(output),
                       const_cast<Real *>(output + output_size),
                       const_cast<Real *>(beginning_valid_window),
                       output,
                       [this](Real const & a, Real const & b) -> Real {
                            return a + b / conv_size;
                        });
        // cblas_daxpy(output_size,             // N    - Number of elements in input vectors
        //             1.0/conv_size,           // DA   - specifies the scalar alpha.
        //             beginning_valid_window,  // DX   - input array
        //             1,                       // INCX - increment for the elements of DX
        //             output,                  // DY,  - output array
        //             1);                      // INCY - increment for the elements of DY

        return *this;
    }

    // void get_full_output(std::vector<Real> *v) const {
    //     if (v->size() != conv_size) {
    //         v->resize(conv_size);
    //     }

    //     std::copy(out, out + conv_size, v->begin());
    // }
};

}}  // end namespace b118::signal
