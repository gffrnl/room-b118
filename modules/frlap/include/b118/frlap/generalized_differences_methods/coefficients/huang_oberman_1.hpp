/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/huang_oberman_1.hpp
 *
 *   Huang & Oberman coefficients - linear interpolation
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
struct huang_oberman_1 : generator<Real, huang_oberman_1> {
    using generator<Real, huang_oberman_1>::coeffs;

    explicit huang_oberman_1(std::size_t n)
        : generator<Real, huang_oberman_1>(n)
    {}

    void generate(Real ealpha, Real deltax) {
        constexpr Real       inv_sqrtpi = b118::numbers::inv_sqrtpi_v<Real>;
        std::size_t    const n         = coeffs.size();

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

            coeffs[1] = is_ealpha_one ? ca * ch * (static_cast<Real>(2)
                                             - log(static_cast<Real>(2)))
                                      : ca * ch * coeff_k_1(ealpha);

            for (std::size_t k = 2; k < n; ++k) {
                coeffs[k] = ca * pow(deltax * static_cast<Real>(k), -ealpha)
                               * series(k, ealpha) / static_cast<Real>(k);
            }
        }
    }

 public:
    static Real series(std::size_t k, Real ealpha) {
        Real const inv_k = 1 / static_cast<Real>(k);
        Real const x     = inv_k * inv_k;
        Real       prod  = 1;
        Real       sum   = 0;
        for (unsigned n = 4; n < 300; n += 2) {
            // prod *= x*(ealpha - 3 + n)*(ealpha - 2 + n)/(n * (n - 1));
            prod *= ((ealpha - 3)/n + 1)*((ealpha - 2)/n + 1)
                    / ((1 - static_cast<Real>(1)/n)) * x;
            auto old_sum = sum;
            sum = old_sum + prod;
            if (sum == old_sum) break;
        }
        return 1 + sum;
    }  // Asymptotics = 1 as k -> infinity

    // The following function calculates the first coefficient
    // for alpha != 1. Though it is not defined for alpha = 1,
    // its limit as alpha -> 1 is 2 - log(2).
    static Real coeff_k_1(Real alpha) {
        Real const l2  = log(static_cast<Real>(2));
        Real const T_1 = expm1(l2*(1 - alpha))/(alpha - 1);
        return 1/(2 - alpha) + (1 + T_1)/alpha;
    }
};

}  // end namespace coefficients
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
