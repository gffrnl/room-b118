/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/gorenflo_mainardi.hpp
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
#include "./generator.hpp"

namespace b118         {
namespace frlap        {
namespace gdm          {
namespace coefficients {

template<typename Real>
struct gorenflo_mainardi : generator<Real, gorenflo_mainardi> {
    using generator<Real, gorenflo_mainardi>::coeffs;

    void generate(Real ealpha, Real deltax) {
        constexpr Real       pi = b118::numbers::pi_v<Real>;
        std::size_t    const n  = coeffs.size();

        //   - there is a way to mitigate the large errors
        //     when alpha ~~ 1 but alpha != 1 ???
        //   - when alpha == 1 it is better use 1.0/(k*(k+1))
        //     or take the logs and then exponentiate ???

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
            if (b118::almost_equal<Real>(ealpha, 1)) {
                Real c1 = static_cast<Real>(1) / (deltax * pi);
                coeffs[0] = c1 * static_cast<Real>(2);

                for (std::size_t k = 1; k < n; ++k)
                    coeffs[k] = - c1
                        / (static_cast<Real>(k) * static_cast<Real>(k + 1));
                // or it is better to use ??
                // coeffs[k] = - c1 / exp(  log(static_cast<Real>(k))
                //             + log(static_cast<Real>(k+1)) );
                //
            } else {
                Real ch = static_cast<Real>(1) / pow(deltax, ealpha);
                Real cc = static_cast<Real>(1)
                    / cos(ealpha * pi / static_cast<Real>(2));
                Real ca = ealpha;

                // REMARK: Using Log-Gamma function and then
                //         exponentiate seems to be better to
                //         mitigate round-off/approximation errors.
                //         We only can do this because the result
                //         of our Gammas are positive, since lgamma()
                //         computes the logarithm of Gamma's ABSOLUTE
                //         value.

                if (ealpha < static_cast<Real>(1)) {
                    coeffs[0] = cc * ch;

                    ca *= (ealpha - static_cast<Real>(1)) * ch
                                / static_cast<Real>(2);
                    for (std::size_t k = 1; k < n; ++k)
                    coeffs[k] = cc * (ca * exp(lgamma(static_cast<Real>(k)
                                                    - ealpha)
                                             - lgamma(static_cast<Real>(k + 1))
                                             - lgamma(static_cast<Real>(2)
                                                    - ealpha)));
                } else {  // 1 < ealpha < 2
                    coeffs[0] = - cc * (ca * ch);

                    ca *= (ealpha - static_cast<Real>(1))
                                        / static_cast<Real>(2);
                    coeffs[1] = cc * ( (static_cast<Real>(1) + ca) * ch)
                                        / static_cast<Real>(2);

                    ca *= (static_cast<Real>(2) - ealpha) * ch;
                    for (std::size_t k = 2; k < n; ++k)
                        coeffs[k] = cc * (ca * exp(lgamma(
                                                       static_cast<Real>(k + 1)
                                                     - ealpha)
                                                 - lgamma(
                                                       static_cast<Real>(k + 2))
                                                 - lgamma(
                                                       static_cast<Real>(3)
                                                     - ealpha)));
                }
            }
        }
    }
};

}  // end namespace coefficients
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
