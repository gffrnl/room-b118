/*   libb118
 *
 *   modules/base/b118/fftw/alloc.hpp
 *
 *   Wrappers for fftw3 fftw_malloc
 *
 *   Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>
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
#include "./complex.hpp"
#include "./complex_ptr.hpp"

namespace b118 {
namespace fftw {
    template<typename Real>
    Real * alloc(std::size_t n) {
        return static_cast<Real *>(fftwq_malloc(sizeof(Real) * n));
    }

    template<>
    float * alloc<float>(std::size_t n) {
        return static_cast<float *>(fftwf_malloc(sizeof(float) * n));
    }

    template<>
    double * alloc<double>(std::size_t n) {
        return static_cast<double *>(fftw_malloc(sizeof(double) * n));
    }

    template<>
    long double * alloc<long double>(std::size_t n) {
        return static_cast<long double *>(fftwl_malloc(sizeof(long double)*n));
    }

    // template<typename Real>
    // typename fftw::complex_ptr<Real>::raw_type
    //     alloc_complex(std::size_t n) {
    //         return static_cast<b118::fftw::complex_ptr<double>::raw_type>(
    //             fftwq_malloc(sizeof(b118::fftw::complex<double>::raw_type)* n));
    // }

    // template<>
    // typename fftw::complex_ptr<float>::raw_type
    //     alloc_complex<float>(std::size_t n) {
    //         return static_cast<b118::fftw::complex_ptr<float>::raw_type>(
    //             fftwf_malloc(sizeof(b118::fftw::complex<float>::raw_type) * n));
    // }

    // template<>
    // typename fftw::complex_ptr<double>::raw_type
    //     alloc_complex<double>(std::size_t n) {
    //         return static_cast<b118::fftw::complex_ptr<double>::raw_type>(
    //             fftw_malloc(sizeof(b118::fftw::complex<double>::raw_type) * n));
    // }

    // template<>
    // typename fftw::complex_ptr<long double>::raw_type
    //     alloc_complex<long double>(std::size_t n) {
    //         return static_cast<b118::fftw::complex_ptr<long double>::raw_type>(
    //             fftw_malloc(sizeof(b118::fftw::complex<long double>::raw_type)
    //                 * n));
    // }

    template<typename Real>
    fftw::complex_ptr<Real> alloc_complex(std::size_t n) {
        fftw::complex_ptr<Real> ptr;
        ptr.raw =
            static_cast<b118::fftw::complex_ptr<double>::raw_type>(
                fftwq_malloc(sizeof(b118::fftw::complex<double>::raw_type)
                             * n));
        return ptr;
    }

    template<>
    fftw::complex_ptr<float> alloc_complex<float>(std::size_t n) {
        fftw::complex_ptr<float> ptr;
        ptr.raw =
            static_cast<b118::fftw::complex_ptr<float>::raw_type>(
                fftwf_malloc(sizeof(b118::fftw::complex<float>::raw_type)
                             * n));
        return ptr;
    }

    template<>
    fftw::complex_ptr<double> alloc_complex<double>(std::size_t n) {
        fftw::complex_ptr<double> ptr;
        ptr.raw =
            static_cast<b118::fftw::complex_ptr<double>::raw_type>(
                fftw_malloc(sizeof(b118::fftw::complex<double>::raw_type)
                            * n));
        return ptr;
    }

    template<>
    fftw::complex_ptr<long double> alloc_complex<long double>(std::size_t n) {
        fftw::complex_ptr<long double> ptr;
        ptr.raw =
            static_cast<b118::fftw::complex_ptr<long double>::raw_type>(
                fftwl_malloc(sizeof(b118::fftw::complex<long double>::raw_type)
                             * n));
        return ptr;
    }

}  // end namespace fftw
}  // end namespace b118
