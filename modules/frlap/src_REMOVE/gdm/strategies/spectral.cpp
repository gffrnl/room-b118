/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/spectral.cpp
 *
 *   Spectral strategy
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
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <b118/frlap/gdm/strategies/spectral.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

#include <iostream>

using b118::frlap::gdm::strategies::spectral;

void spectral::generate_coefficients(double      ealpha,
                                     double      deltax,
                                     double*     coeffs,
                                     std::size_t n) const {
  using boost::math::hypergeometric_pFq;

  {  // Treat the boundary cases ealpha = 0 and ealpha = 2
    if (b118::almost_equal(ealpha, 0.0)) {
      coeffs[0] = 1.0;
      std::fill(coeffs + 1, coeffs + n, 0.0);
      return;
    }

    if (b118::almost_equal(ealpha, 2.0)) {
      double const c2 = 2.0 / (deltax * deltax);
      coeffs[0] =  (M_PI * M_PI) / (3.0 * deltax * deltax);
      for (std::size_t k = 1; k < n; k+=2)
        coeffs[k] = - c2 / (k * k);
      for (std::size_t k = 2; k < n; k+=2)
        coeffs[k] =   c2 / (k * k);
      return;
    }
  }

  { // Treat the case ealpha = 1
    if (b118::almost_equal(ealpha, 1.0)) {
      double const c1 = 1.0 / (M_PI * deltax);
      coeffs[0] =  M_PI / (2 * deltax);
      for (std::size_t k = 1; k < n; k+=2)
        coeffs[k] = - c1 * (2.0 / (k * k));
      for (std::size_t k = 2; k < n; k+=2)
        coeffs[k] = 0.0;
      return;
    }
  }

  { // Treat the cases {0 < ealpha < 1} U {1 < ealpha < 2}
    double const a1 = (ealpha + static_cast<double>(1)) / 2;
    double const b1 = static_cast<double>(1) / 2;
    double const b2 = (ealpha + static_cast<double>(3)) / 2;
    double const ch = std::pow(M_PI / deltax, ealpha) / (ealpha + 1);
    double const mult = - M_PI * M_PI / static_cast<double>(4);
    for (std::size_t k = 1; k < n; ++k) {
      coeffs[k] = ch * hypergeometric_pFq({a1}, {b1, b2}, mult * (k * k));
      if (k == 1)
        std::cout << "spec coeffs[" <<k <<"] = " << coeffs[1] << std::endl;
    }
  }
}
