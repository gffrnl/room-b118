/*   libb118
 *
 *   modules/base/b118/almost_equal.hpp
 *     Contains:
 *       - Facilities for comparison of floating-point numbers
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

#include <cstdint>
#include <cmath>
#include <limits>

namespace b118 {

    // Almost equal
    // Implementation of numpy `within_tol`
    // <https://github.com/numpy/numpy/blob/v1.26.0/numpy/core/numeric.py#L2330-L2333>
    // TODO(fazedo): change atol and rtol to fit the type
    template<typename Real>
    inline bool almost_equal(Real x, Real y,
                             Real atol =
                                sqrt(std::numeric_limits<Real>::epsilon()),
                             Real rtol =
                                cbrt(std::numeric_limits<Real>::epsilon())) {
        return std::fabs(x - y) < atol + rtol * std::fabs(y);
    }

    // Almost equal with `dig` base-10 significant digits
    // Reference: R.L. Burden and J.D. Faires. Análise Numérica.
    //            Cengage Learnig, 8 edition, 2013.
    template<typename Real>
    inline bool almost_equal_dig(Real x, Real y,
                                 std::uint16_t dig =
                                     std::numeric_limits<Real>::digits10 - 1) {
        return std::fabs(x - y) <
            (static_cast<Real>(5) /
                std::pow(static_cast<Real>(10) ,
                         static_cast<Real>(dig)) * std::abs(y));
    }

}  // namespace b118
