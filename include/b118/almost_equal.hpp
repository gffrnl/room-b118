/*   libb118
 *
 *   include/b118/almost_equal.hpp
 *     Contains:
 *       - Facilities for comparison of floating-point numbers
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

#include <cmath>
#include <functional>

namespace b118 {


    // Almost equal (Python within_tol)
    // TODO: change atol and rtol to fit the type
    template<typename Real>
    inline bool almost_equal(Real x, Real y, Real atol = 1.0e-5, Real rtol = 1.0e-8) {
        return std::less_equal{}(std::fabs(x - y), atol + rtol * std::fabs(y));
    }
    
    // Almost equal with `dig` base-10 significant digits
    // Reference: R.L. Burden and J.D. Faires. Análise Numérica. Cengage Learning, 8 edition, 2013.
    template<typename Real>
    inline bool almost_equal_dig(Real x, Real y, short dig = std::numeric_limits<T>::digits10 - 1) {
        return std::less{}(
            std::fabs(x / y - 1),
            static_cast<Real>(5) / std::pow(static_cast<Real>(10), static_cast<Real>(dig))
        );
    }
}