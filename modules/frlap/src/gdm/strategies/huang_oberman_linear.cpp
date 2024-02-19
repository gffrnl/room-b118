/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/huang_oberman_linear.cpp
 *   
 *   Huang & Oberman linear strategy
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

#include <cstddef>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <b118/almost_equal.hpp>
#include <b118/numbers.hpp>
#include <b118/frlap.hpp>
#include <b118/frlap/gdm/strategies/huang_oberman_linear.hpp>

extern double d1G_alpha_ne_1(double, std::size_t);
extern double d2G_alpha_ne_1(double, std::size_t);
extern double d1G_alpha_eq_1(std::size_t);
extern double d2G_alpha_eq_1(std::size_t);

using b118::frlap::gdm::strategies::huang_oberman_linear;

void huang_oberman_linear::generate_coefficients(double      ealpha,
                                                 double      deltax,
                                                 double*     coeffs,
                                                 std::size_t n) const {
  double constexpr inv_sqrtpi = b118::numbers::inv_sqrtpi_v<double>;

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
    bool const is_ealpha_one = b118::almost_equal(ealpha, 1.0);

    double const ca = is_ealpha_one ? - frlap::normal_const<1>(1.0)
                                    : - frlap::normal_const<1>(ealpha);

    double const ch = is_ealpha_one ? + 1.0 / deltax
                                    : + std::pow(deltax, -ealpha);

    // REMARK: using Log-Gamma function and then
    //         exponentiate seems better to mitigate
    //         round-off/approximation errors.
    //         Since lgamma() computes the log of
    //         the abs value of Gamma, we only can
    //         do this because the arguments of
    //         our Gammas are always positive.
    //
    coeffs[0] = ch * inv_sqrtpi * std::exp2(ealpha)
                   * std::exp(std::lgamma(0.5 + 0.5 * ealpha)
                            - std::lgamma(2.0 - 0.5 * ealpha));

    double const ca_ch = ca * ch;

    if (is_ealpha_one == true) {  // Could test only ealpha == 1
      coeffs[1] =   ca_ch * (1.0 - d2G_alpha_eq_1(1)
                                 + d1G_alpha_eq_1(2)
                                 - d1G_alpha_eq_1(1));

      for (std::size_t k = 2; k < n; ++k)
        coeffs[k] = ca_ch * (        d1G_alpha_eq_1(k + 1)  // NOLINT
                             - 2.0 * d1G_alpha_eq_1(k)
                             +       d1G_alpha_eq_1(k - 1));
    } else {
      coeffs[1] =   ca_ch * (1.0 / (2.0 - ealpha) - d2G_alpha_ne_1(ealpha, 1)
                                                  + d1G_alpha_ne_1(ealpha, 2)
                                                  - d1G_alpha_ne_1(ealpha, 1));

      for (std::size_t k = 2; k < n; ++k)
        coeffs[k] = ca_ch * (        d1G_alpha_ne_1(ealpha, k+1)  // NOLINT
                             - 2.0 * d1G_alpha_ne_1(ealpha, k)
                             +       d1G_alpha_ne_1(ealpha, k-1));
    }
  }
}
