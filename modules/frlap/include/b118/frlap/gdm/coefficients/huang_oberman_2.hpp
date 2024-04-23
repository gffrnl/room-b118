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
#include "./method.hpp"
#include "./generator.hpp"

namespace b118 { namespace frlap { namespace gdm { namespace coefficients {

template<typename Real>
class huang_oberman_2 : public method<Real>,
                        public generator<Real, huang_oberman_2> {
    using generator<Real, huang_oberman_2>::ealpha;
    using generator<Real, huang_oberman_2>::deltax;

 public:
    huang_oberman_2(Real const & ealpha, Real const & deltax)
        : generator<Real, huang_oberman_2>(ealpha, deltax),
          is_ealpha_one(b118::almost_equal<Real>(ealpha, 1)),
          ca(is_ealpha_one ? - frlap::normal_const<1, Real>(1)
                           : - frlap::normal_const<1, Real>(ealpha)),
          ch(is_ealpha_one ? + static_cast<Real>(1) / deltax
                           : + pow(deltax, -ealpha))
    {}

    Real operator()(std::size_t const & k) const override {
        // The boundary cases ealpha = 0 and ealpha = 2
        if (b118::almost_equal<Real>(ealpha, 0)) {
            if (k == 0) return static_cast<Real>(1);
            return static_cast<Real>(0);
        }
        if (b118::almost_equal<Real>(ealpha, 2)) {
            if (k == 0)
                return   static_cast<Real>(2) / (deltax * deltax);
            if (k == 1)
                return - static_cast<Real>(1) / (deltax * deltax);
            return static_cast<Real>(0);
        }

        // The general case 0 < ealpha < 2
        // REMARK: using Log-Gamma function and then
        //         exponentiate seems better to mitigate
        //         round-off/approximation errors.
        //         Since lgamma() computes the log of
        //         the abs value of Gamma, we only can
        //         do this because the arguments of
        //         our Gammas are always positive.
        //
        if (k == 0)
            return ch * b118::numbers::inv_sqrtpi_v<Real> * exp2(ealpha)
                           * exp(lgamma((static_cast<Real>(1) + ealpha)
                                               / static_cast<Real>(2))
                               - lgamma(static_cast<Real>(2)
                                      - ealpha / static_cast<Real>(2)));

        if (k == 1) {
            // Remark: coeffs[1] is continuous at ealpha = 1 but we prefered to
            //         treat it separately.
            if (is_ealpha_one == true)
                return ca * ch * (static_cast<Real>(8)
                                     - static_cast<Real>(5)
                                            * log(static_cast<Real>(3)))
                                      / static_cast<Real>(2);

            return ca * ch * coeff_k_1(ealpha);
        }

        // The formulation in Taylor series allows us to use
        // the same loop for both cases: alpha = 1 e alpha != 1
        if (k%2 == 0)  // k even
            return 2.0 * ca * series_even(k, ealpha)
                                     * pow(deltax * static_cast<Real>(k),
                                         - ealpha)
                                     / static_cast<Real>(k);

        // k odd
        return ca * series_odd(k, ealpha)
                  * pow(deltax * static_cast<Real>(k),
                        - ealpha)
                  / static_cast<Real>(k);
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

 private:
    bool const is_ealpha_one;
    Real const ca;
    Real const ch;
};

}}}}  // end namespace b118::frlap::gdm::coefficients
