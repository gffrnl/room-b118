/*   libb118
 *
 *   modules/base/b118/fftw/plan.hpp
 *
 *   Wrappers for fftw3 plans
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
struct plan final {
    fftwq_plan raw;
    plan() : raw(nullptr) {}
};

template<>
struct plan<float> final {
    fftwf_plan raw;
    plan() : raw(nullptr) {}
};

template<>
struct plan<double> final {
    fftw_plan raw;
    plan() : raw(nullptr) {}
};

template<>
struct plan<long double> final {
    fftwl_plan raw;
    plan() : raw(nullptr) {}
};

template<typename Real>
void destroy_plan(plan<Real> & p) {  // NOLINT
    fftwq_destroy_plan(p.raw);
}

template<>
void destroy_plan<float>(plan<float> & p) {  // NOLINT
    fftwf_destroy_plan(p.raw);
}

template<>
void destroy_plan<double>(plan<double> & p) {  // NOLINT
    fftw_destroy_plan(p.raw);
}

template<>
void destroy_plan<long double>(plan<long double> & p) {  // NOLINT
    fftwl_destroy_plan(p.raw);
}


}  // end namespace fftw
}  // end namespace b118
