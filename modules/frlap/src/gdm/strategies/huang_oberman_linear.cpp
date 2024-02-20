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


template <typename Real>
Real series(int k, Real alpha) {
  Real inv_k = 1/static_cast<Real>(k);
  Real x = inv_k*inv_k;
  Real prod = 1/static_cast<Real>(2);
  Real sum = 0;
  for (int n = 4 ; n < 300; n += 2) {
      prod *= x*(alpha - 3 + n)*(alpha - 2 + n)/(n * (n - 1));
          auto old_sum = sum;
      sum = old_sum + prod;
      if (sum == old_sum) {
          break;
      }
  }
  return (1 + 2*sum);
}  // Asymptotics = 1 as k -> infinity


// The following function calculates the first coefficient
// for alpha != 1. Though it is not defined for alpha = 1,
// its limit as alpha -> 1 is 2 - log(2).
template <typename Real>
Real coeff_k_1(Real alpha) {
  Real l2 = std::log(static_cast<Real>(2));
  Real T_1 = std::expm1(l2*(1 - alpha))/(alpha - 1);
  return 1/(2 - alpha) + (1 + T_1)/alpha;
}


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

    coeffs[1] = is_ealpha_one ? ca * ch * (2.0 - std::log(2.0))
                              : ca * ch * coeff_k_1(ealpha);

    for (std::size_t k = 2; k < n; ++k) {
      coeffs[k] = ca * std::pow(k*deltax, -ealpha) * series(k, ealpha) / k;
    }
  }
}
