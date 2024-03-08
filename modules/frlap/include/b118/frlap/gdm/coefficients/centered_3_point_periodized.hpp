/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/centered_3_point_periodized.hpp
 *
 *   Coefficients by periodization of 3-point rule coefficients
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

namespace b118         {
namespace frlap        {
namespace gdm          {
namespace coefficients {

template<typename Real>
struct centered_3_point_periodized :
    generator<Real, centered_3_point_periodized> {
    using generator<Real, centered_3_point_periodized>::coeffs;

    void generate(Real ealpha, Real deltax) {
        constexpr Real       pi = b118::numbers::pi_v<Real>;
        std::size_t    const n  = coeffs.size();

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
            Real const ch = static_cast<Real>(1)/ (pi * pow(deltax, ealpha));
            Real const cs = - ch * sin(ealpha * pi / static_cast<Real>(2));

            coeffs[0] = ch * (exp2(ealpha)
                * std::beta((ealpha + 1) / static_cast<Real>(2),
                            static_cast<Real>(1) / static_cast<Real>(2)));

            for (std::size_t k = 1; k < n; ++k)
            coeffs[k] = cs * std::beta(static_cast<Real>(k)
                                     - ealpha / static_cast<Real>(2),
                                       ealpha + static_cast<Real>(1));
        }
    }
};

}  // end namespace coefficients
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
