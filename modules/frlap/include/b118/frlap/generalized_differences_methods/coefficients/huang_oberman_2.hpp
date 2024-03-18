/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/huang_oberman_2.hpp
 *
 *   Huang & Oberman coefficients - quadratic interpolation
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

#include <cstddef>
#include <vector>
#include <b118/almost_equal.hpp>
#include <b118/numbers.hpp>
#include "../../../frlap.hpp"
#include "./generator.hpp"

namespace b118         {
namespace frlap        {
namespace gdm          {
namespace coefficients {

template<typename Real>
struct huang_oberman_2 : generator<Real, huang_oberman_2> {
    using generator<Real, huang_oberman_2>::coeffs;

    explicit huang_oberman_2(std::size_t n)
        : generator<Real, huang_oberman_2>(n)
    {}

    void generate(Real ealpha, Real deltax) {
        constexpr Real       inv_sqrtpi = b118::numbers::inv_sqrtpi_v<Real>;
        std::size_t    const n          = coeffs.size();

        {  // Treat the boundary cases ealpha = 0 and ealpha = 2
            if (b118::almost_equal<Real>(ealpha, 0)) {
                coeffs[0] = static_cast<Real>(1);
                std::fill(coeffs.begin() + 1, coeffs.end(),
                          static_cast<Real>(0));
                return;
            }

            if (b118::almost_equal<Real>(ealpha, 2)) {
                coeffs[0] =   static_cast<Real>(2) / (deltax * deltax);
                coeffs[1] = - static_cast<Real>(1) / (deltax * deltax);
                std::fill(coeffs.begin() + 2, coeffs.end(),
                          static_cast<Real>(0));
                return;
            }
        }

        {  // Now we treat the general case 0 < ealpha < 2
            bool const is_ealpha_one
                          = b118::almost_equal<Real>(ealpha, 1);

            Real const ca = is_ealpha_one
                          ? - frlap::normal_const<1, Real>(1)
                          : - frlap::normal_const<1, Real>(ealpha);

            Real const ch = is_ealpha_one
                          ? + static_cast<Real>(1) / deltax
                          : + pow(deltax, -ealpha);

            // REMARK: using Log-Gamma function and then
            //         exponentiate seems better to mitigate
            //         round-off/approximation errors.
            //         Since lgamma() computes the log of
            //         the abs value of Gamma, we only can
            //         do this because the arguments of
            //         our Gammas are always positive.
            //
            coeffs[0] = ch * inv_sqrtpi * exp2(ealpha)
                           * exp(lgamma((static_cast<Real>(1) + ealpha)
                                               / static_cast<Real>(2))
                               - lgamma(static_cast<Real>(2)
                                      - ealpha / static_cast<Real>(2)));

            // Remark: coeffs[1] is continuous at ealpha = 1 but we prefered to
            //         treat it separately.

            if (is_ealpha_one == true) {  // Could test only ealpha == 1
                coeffs[1] = ca * ch * (static_cast<Real>(8)
                                     - static_cast<Real>(5)
                                            * log(static_cast<Real>(3)))
                                      / static_cast<Real>(2);
            } else {
                coeffs[1] = ca * ch * coeff_k_1(ealpha);
            }

            // The formulation in Taylor series allows us to use
            // the same loop for both cases: alpha = 1 e alpha != 1
            for (std::size_t k = 2; k < n; k += 2)
                coeffs[k] = 2.0 * ca * series_even(k, ealpha)
                                     * pow(deltax * static_cast<Real>(k),
                                         - ealpha)
                                     / static_cast<Real>(k);

            for (std::size_t k = 3; k < n; k += 2)
                coeffs[k] =       ca * series_odd(k, ealpha)
                                     * pow(deltax * static_cast<Real>(k),
                                         - ealpha)
                                     / static_cast<Real>(k);
        }
    }

 public:
    // The following two functions are involved in calculating the coefficients
    // The formulation via Taylor Series is regular at alpha = 1 and removes
    // cathastrophic cancelations near both alpha = 1 and k -> infinity.
    // In both the first term is summed in the end, so we have a bit of
    // incresed accuracy.
    // Both present the asymptotics 4/6 as k -> infinity.
    // This asymptotic could be absorved in the main loop.
    // Obs1: The series are truncated at 300 terms, which is more than enough
    //      for 50 digits.

    static Real series_even(std::size_t k, Real ealpha) {
        Real constexpr assymp  = static_cast<Real>(4) / 6;
        Real const     inv_k   = 1 / static_cast<Real>(k);
        Real const     x       = inv_k * inv_k;
        Real const     T0      = 1;  // first term of the series
        Real           prod    = static_cast<Real>(1) / 2;
        Real           sum     = 0;
        Real           old_sum;

        for (unsigned n = 4; n < 300; n += 2) {
            prod *= ((ealpha - 3)/n + 1) * ((ealpha - 2)/n + 1)
                    / ((1 + static_cast<Real>(1)/n)) * x;
            old_sum = sum;
            sum = old_sum + prod * n;
            if (sum == old_sum) break;
        }
        return (T0 + sum) * assymp;
    }  // Asymptotics = 4/6 as k -> infinity

    static Real series_odd(std::size_t k, Real ealpha) {
        Real constexpr assymp = static_cast<Real>(4) / 6;
        Real const     inv_2k = 2 / static_cast<Real>(k);
        Real const     x      = inv_2k * inv_2k;
        Real const     T0     = 1;
        Real           sum    = T0;
        Real           prod   = T0;
        Real           old_sum;

        for (unsigned n = 4 ; n < 300; n += 2) {
            prod *= ((ealpha - 3)/n + 1)*((ealpha - 2)/n + 1)
                    / ((1 + static_cast<Real>(1)/n)) * x;
            old_sum = sum;
            sum = old_sum - prod * (n - 3);
            if (sum == old_sum) break;
        }
        return sum * assymp;
    }

    // This function calcultes the coefficient para k = 1 for ealpha != 1.
    // The limit for ealpha -> 1 is (8 - 5*std::log(3))/2,
    // which comes from:
    // (exp(a*x) - 1)/x = a,   x -> 0.
    static Real coeff_k_1(Real ealpha) {  // k = 1 alpha != 1
        Real const log3 = log(3);
        Real const T_1  = - expm1(log3*(1 - ealpha)) / (1 - ealpha);
        Real const T_2  = 4 + (ealpha + 4) * T_1 / 2;
        return T_2 / (ealpha * (2 - ealpha));
    }
};

}  // end namespace coefficients
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
