/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/huang_oberman_quadratic.cpp
 *   
 *   Huang & Oberman quadratic strategy
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
#include <b118/frlap/gdm/strategies/huang_oberman_quadratic.hpp>


// TODO(Guilherme): Better names for the functions?
// TODO(Guilherme): Is it necessary to use almost_equal in this file?
//                : All expression are well conditioned at alpha = 1.


// The following two functions are involved in calculating the coefficients
// The formulation via Taylor Series is regular at alpha = 1 and removes
// cathastrophic cancelations near both alpha = 1 and k -> infinity.
// In both the first term is summed in the end, so we have a bit of increased
// accuracy.
// Both present the asymptotics 4/6 as k -> infinity.
// This asymptotic could be absorved in the main loop.
// Obs1: The series are truncated at 300 terms, which is more than enough
//      for 50 digits.
template<class Real>
Real series_even(int k, Real ealpha) {
  Real constexpr assymp = static_cast<Real>(4) / 6;
  Real const inv_k = 1/static_cast<Real>(k);
  Real const x = inv_k*inv_k;
  Real const T0 = 1;  // first term of the series
  Real prod = static_cast<Real>(1)/2;
  Real sum = 0;
  Real old_sum;

  for (int n = 4 ; n < 300; n += 2) {
    prod *= ((ealpha - 3)/n + 1)*((ealpha - 2)/n + 1)
            / ((1 + static_cast<Real>(1)/n)) * x;
    old_sum = sum;
    sum = old_sum + prod * n;
    if (sum == old_sum) {
      break;
    }
  }
  return (T0 + sum) * assymp;
}  // Asymptotics = 4/6 as k -> infinity

template<class Real>
Real series_odd(int k, Real ealpha) {
  Real constexpr assymp = static_cast<Real>(4) / 6;
  Real const inv_2k = 2 / static_cast<Real>(k);
  Real const x = inv_2k*inv_2k;
  Real const T0 = 1;
  Real sum = T0;
  Real prod = T0;
  Real old_sum;

  for (int n = 4 ; n < 300; n += 2) {
      prod *= ((ealpha - 3)/n + 1)*((ealpha - 2)/n + 1)
            / ((1 + static_cast<Real>(1)/n)) * x;
      old_sum = sum;
      sum = old_sum + prod * (3 - n);
      if (sum == old_sum) {
          break;
      }
  }
  return sum * assymp;
}

// This function calcultes the coefficient para k = 1 for ealpha != 1.
// The limit for ealpha -> 1 is (8 - 5*std::log(3))/2,
// which comes from:
// (exp(a*x) - 1)/x = a,   x -> 0.
template<class Real>
Real better_1(Real ealpha) {  // k = 1 alpha != 1
  Real constexpr log3 = std::log(3);
  Real const T_1 = -std::expm1(log3*(1 - ealpha)) / (1 - ealpha);
  Real const T_2 = 4 + (ealpha + 4) * T_1 / 2;
  return T_2 / (ealpha * (2 - ealpha));
}

// extern double d0G_alpha_ne_1(double, std::size_t);
// extern double d1G_alpha_ne_1(double, std::size_t);
// extern double d2G_alpha_ne_1(double, std::size_t);
// extern double d0G_alpha_eq_1(std::size_t);
// extern double d1G_alpha_eq_1(std::size_t);
// extern double d2G_alpha_eq_1(std::size_t);


using b118::frlap::gdm::strategies::huang_oberman_quadratic;

void huang_oberman_quadratic::generate_coefficients(double      ealpha,
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

    // Remark: coeffs[1] is continuous at ealpha = 1 but we prefered to
    //         treat it separately.

    if (is_ealpha_one == true) {  // Could test only ealpha == 1
      coeffs[1] = ca * ch * (8 - 5*std::log(3))/2;
    } else {
      coeffs[1] = ca * ch * better_1(ealpha);
    }

    // The formulation in Taylor series allows us to use
    // the same loop for both cases: alpha = 1 e alpha != 1
    for (std::size_t k = 2; k < n; k += 2)
      coeffs[k] = 2.0 * ca * series_even(k, ealpha)
                           * std::pow(deltax * k, - ealpha) / k;

    for (std::size_t k = 3; k < n; k += 2)
      coeffs[k] =       ca * series_odd(k, ealpha)
                           * std::pow(deltax * k, - ealpha) / k;
  }
}
