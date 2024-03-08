/*   libb118
 *
 *   modules/base/b118/fftw/complex.hpp
 *
 *   Wrappers for fftw3 complex
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

namespace b118 {
namespace fftw {

template<typename Real>
struct complex final {
    using raw_type = fftwq_complex;
    raw_type raw;
    complex() : raw({static_cast<Real>(0), static_cast<Real>(0)}) {}
};

template<>
struct complex<float> final {
    using raw_type = fftwf_complex;
    raw_type raw;
    complex() : raw({0.0f, 0.0f}) {}
};

template<>
struct complex<double> final {
    using raw_type = fftw_complex;
    raw_type raw;
    complex() : raw({0.0, 0.0}) {}
};

template<>
struct complex<long double> final {
    using raw_type = fftwl_complex;
    raw_type raw;
    complex() : raw({0.0l, 0.0l}) {}
};


}  // end namespace fftw
}  // end namespace b118
