/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/huang_oberman_quadratic.cpp
 *   
 *   Huang & Oberman Quadratic strategy
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
#include <cfloat>
#include <b118/almost_equal.hpp>
#include <b118/frlap.hpp>
#include <b118/frlap/gdm/strategies/huang_oberman_quadratic.hpp>

#ifndef M_1_SQRTPI
#define M_1_SQRTPI 0.564189583547756286948079451561
#endif

extern double d0G_alpha_ne_1(double, std::size_t);
extern double d1G_alpha_ne_1(double, std::size_t);
extern double d2G_alpha_ne_1(double, std::size_t);
extern double d0G_alpha_eq_1(std::size_t);
extern double d1G_alpha_eq_1(std::size_t);
extern double d2G_alpha_eq_1(std::size_t);


using b118::frlap::gdm::strategies::huang_oberman_quadratic;

void huang_oberman_quadratic::generate_coefficients(double      ealpha,
                                                    double      deltax,
                                                    double*     coeffs,
                                                    std::size_t n) const {
    {  // treat the boundary cases ealpha = 0 and ealpha = 2
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

    {  // {0 < ealpha < 2}
        double ca = 0.0;
        double ch = 0.0;

        if (b118::almost_equal(ealpha, 1.0)) {
            ca = - frlap::normal_const<1>(1.0);
            ch = 1.0 / deltax;
        } else {
            ca = - frlap::normal_const<1>(ealpha);
            ch = std::pow(deltax, -ealpha);
        }


        // REMARK: using Log-Gamma function and then
        //         exponentiate seems better to mitigate
        //         round-off/approximation errors.
        //         Since lgamma() computes the log of
        //         the abs value of Gamma, we only can
        //         do this because the arguments of
        //         our Gammas are always positive.
        //
        coeffs[0] = ch * M_1_SQRTPI * std::exp2(ealpha) *
            std::exp(std::lgamma(0.5 * (ealpha + 1.0))
                - std::lgamma(2.0 -  0.5 * ealpha));

        ch *= ca;
        if (b118::almost_equal(ealpha, 1.0)) {
            coeffs[1] = ch * (1.0
                -       d2G_alpha_eq_1(1)
                - 0.5 * d1G_alpha_eq_1(3)
                - 1.5 * d1G_alpha_eq_1(1)
                +       d0G_alpha_eq_1(3)
                -       d0G_alpha_eq_1(1));

            for (std::size_t k = 2; k < n; k+=2)
                coeffs[k] = ch * 2.0 * (d1G_alpha_eq_1(k+1)
                        + d1G_alpha_eq_1(k-1)
                                        - d0G_alpha_eq_1(k+1)
                                        + d0G_alpha_eq_1(k-1));

            for (std::size_t k = 3; k < n; k+=2)
                coeffs[k] = ch * (- 0.5 * (d1G_alpha_eq_1(k+2)
                            + d1G_alpha_eq_1(k-2))
                        - 3.0 * d1G_alpha_eq_1(k)
                        +       d0G_alpha_eq_1(k+2)
                        -       d0G_alpha_eq_1(k-2));
        } else {
            coeffs[1] = ch * (1.0 / (2.0 - ealpha)
                        -       d2G_alpha_ne_1(ealpha, 1)
                        - 0.5 * d1G_alpha_ne_1(ealpha, 3)
                        - 1.5 * d1G_alpha_ne_1(ealpha, 1)
                        +       d0G_alpha_ne_1(ealpha, 3)
                        -       d0G_alpha_ne_1(ealpha, 1));

            for (std::size_t k = 2; k < n; k+=2)
                coeffs[k] = ch * 2.0 * (d1G_alpha_ne_1(ealpha, k+1)
                        + d1G_alpha_ne_1(ealpha, k-1)
                        - d0G_alpha_ne_1(ealpha, k+1)
                        + d0G_alpha_ne_1(ealpha, k-1));

            for (std::size_t k = 3; k < n; k+=2)
                coeffs[k] = ch * (- 0.5 * (d1G_alpha_ne_1(ealpha, k+2)
                                        + d1G_alpha_ne_1(ealpha, k-2) )
                                - 3.0 * d1G_alpha_ne_1(ealpha, k)
                                +       d0G_alpha_ne_1(ealpha, k+2)
                                -       d0G_alpha_ne_1(ealpha, k-2));
        }
    }
}
