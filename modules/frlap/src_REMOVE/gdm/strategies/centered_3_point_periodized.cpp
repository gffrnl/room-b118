/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/centered_3_point_periodized.cpp
 *   
 *   3-point centered periodized strategy
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

#include <cstddef>
#include <cmath>
#include <algorithm>
#include <b118/almost_equal.hpp>
#include <b118/frlap/gdm/strategies/centered_3_point_periodized.hpp>

using b118::frlap::gdm::strategies::centered_3_point_periodized;

void centered_3_point_periodized::generate_coefficients(double      ealpha,
                                                        double      deltax,
                                                        double*     coeffs,
                                                        std::size_t n) const {
  {  // Treat the boundary cases ealpha = 0 and ealpha = 2
    if (b118::almost_equal(ealpha, 0.0)) {
        coeffs[0] = 1.0;
        std::fill(coeffs + 1, coeffs + n, 0.0);
        return;
    }

    if (b118::almost_equal(ealpha, 2.0)) {
        coeffs[0] =  2.0 / (deltax * deltax);
        coeffs[1] = -1.0 / (deltax * deltax);
        std::fill(coeffs + 2, coeffs + n, 0.0);
        return;
    }
  }

  {  // Now we treat the general case 0 < ealpha < 2
    double const ch = 1.0 / (M_PI * std::pow(deltax, ealpha));
    double const cs = - ch * std::sin(0.5 * ealpha * M_PI);

    coeffs[0] = ch * (std::exp2(ealpha) * std::beta(0.5 + 0.5 * ealpha, 0.5));

    for (std::size_t k = 1; k < n; ++k)
      coeffs[k] = cs * std::beta(static_cast<double>(k) - 0.5 * ealpha,
                                 1.0 + ealpha);
  }
}
