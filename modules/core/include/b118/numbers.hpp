/*   libb118
 *
 *   include/b118/numbers.hpp
 *   
 *   Mathematical constants following std C++20 implementation <numbers>
 *
 *   Copyright (C) 2024   Guilherme F. Fornel        <gffrnl@gmail.com>
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

#include <b118/type_traits.hpp>

namespace b118 {
namespace numbers {

    // Reference:
    // <https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/include/std/numbers>


    // pi
    template<typename Real>
    inline constexpr Real pi_v
      = b118::enable_if_real_t<Real>(3.141592653589793238462643383279502884L);

    // 1/sqrt(pi)
    template<typename Real>
    inline constexpr Real inv_sqrtpi_v
      = b118::enable_if_real_t<Real>(0.564189583547756286948079451560772586L);

    // for double
    inline constexpr double pi = pi_v<double>;
    inline constexpr double inv_sqrt_pi = inv_sqrtpi_v<double>;

}  // namespace numbers
}  // namespace b118
