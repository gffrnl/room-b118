/*   libb118
 *
 *   modules/base/b118/fftw/alloc.hpp
 *
 *   Wrappers for fftw3 fftw_free
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
    void free(Real * ptr) {
        fftwq_free(ptr);
    }

    template<>
    void free<float>(float * ptr) {
        fftwf_free(ptr);
    }

    template<>
    void free<double>(double * ptr) {
        fftw_free(ptr);
    }

    template<>
    void free<long double>(long double * ptr) {
        fftwl_free(ptr);
    }

    template<typename Real>
    void free_complex(fftw::complex_ptr<Real> & ptr) {                 // NOLINT
        fftwq_free(ptr.raw);
    }

    template<>
    void free_complex(fftw::complex_ptr<float> & ptr) {                // NOLINT
        fftwf_free(ptr.raw);
    }

    template<>
    void free_complex(fftw::complex_ptr<double> & ptr) {               // NOLINT
        fftw_free(ptr.raw);
    }

    template<>
    void free_complex(fftw::complex_ptr<long double> & ptr) {          // NOLINT
        fftwl_free(ptr.raw);
    }
}  // end namespace fftw
}  // end namespace b118
