/*   libb118
 *
 *   modules/base/b118/fftw/complex.hpp
 *
 *   Wrappers for fftw3 complex pointer
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
#include <cstddef>

namespace b118 {
namespace fftw {

template<typename Real>
struct complex_ptr final {
    using raw_type = fftwq_complex *;
    raw_type raw;
    complex_ptr() : raw(nullptr) {}
};

template<>
struct complex_ptr<float> final {
    using raw_type = fftwf_complex *;
    raw_type raw;
    complex_ptr() : raw(nullptr) {}
};

template<>
struct complex_ptr<double> final {
    using raw_type = fftw_complex *;
    raw_type raw;
    complex_ptr() : raw(nullptr) {}
};

template<>
struct complex_ptr<long double> final {
    using raw_type = fftwl_complex *;
    raw_type raw;
    complex_ptr() : raw(nullptr) {}
};

// template<typename Real>
// struct complex_unique_ptr final {
//     using raw_type = fftwq_complex *;
//     raw_type raw;
//     complex_unique_ptr() : raw(nullptr) {}
//     explicit complex_unique_ptr(std::size_t n)
//         : raw(static_cast<fftwq_complex *>(fftwq_malloc(sizeof(fftwq_complex) * n)))
//     {}
//     ~complex_unique_ptr() {
//         fftwq_free(raw);
//     }
// };

// template<>
// struct complex_unique_ptr<double> final {
//     using raw_type = fftw_complex *;
//     raw_type raw;
//     complex_unique_ptr() : raw(nullptr) {}
//     explicit complex_unique_ptr(std::size_t n)
//         : raw(static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * n)))
//     {}
//     ~complex_unique_ptr() {
//         fftw_free(raw);
//     }
// };

}  // end namespace fftw
}  // end namespace b118
