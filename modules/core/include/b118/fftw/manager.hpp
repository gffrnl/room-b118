/*   libb118
 *
 *   modules/base/b118/fftw/manager.hpp
 *
 *   Manager for fftw3
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

#ifdef B118_FFTW_MANAGER_LOG
#include <iostream>
#endif

namespace b118 {
namespace fftw {

template<typename Real>
struct manager final {
    static void cleanup() {
        fftwq_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
        std::clog << "*** fftwq_cleanup() called ***" << std::endl;
#endif
    }
};

template<>
struct manager<float> final {
    static void cleanup() {
        fftwf_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
        std::clog << "*** fftwf_cleanup() called ***" << std::endl;
#endif
    }
};

template<>
struct manager<double> final {
    static void cleanup() {
        fftw_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
        std::clog << "*** fftw_cleanup() called ***" << std::endl;
#endif
    }
};

template<>
struct manager<long double> final {
    static void cleanup() {
        fftwl_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
        std::clog << "*** fftwl_cleanup() called ***" << std::endl;
#endif
    }
};

}  // end namespace fftw
}  // end namespace b118
