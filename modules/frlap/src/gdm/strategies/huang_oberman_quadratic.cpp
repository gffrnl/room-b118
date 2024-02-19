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

// The following two functions are involved in calculating the coefficients
// The formulation via Taylor Series is regular at alpha = 1 and removes
// cathastrophic cancelations near both alpha = 1 and k -> infinity.
// In both the first term is summed in the end, so we have a bit of increased accuracy.
// Both present the asymptotics 4/6 as k -> infinity
template <class Real>
Real serie_par(int k, Real alpha) {
    Real inv_k = 1/static_cast<Real>(k);
    Real x = inv_k*inv_k;
    Real x_2n = 1;
    Real prod = 1/static_cast<Real>(6);  // (alpha-1)*alpha/(n+1)! para n=2.
    Real T0 = 2*prod;  // first term of the series
    Real soma = 0;

    for (int n = 4 ; n < 300; n += 2) {
        prod *= (alpha - 3 + n)*(alpha - 2 + n)/(n*(n+1));
        x_2n *= x;
        auto soma_old = soma;
        soma = soma_old + prod*x_2n*n;
        if (soma == soma_old) {
            break;
        }
    }
    return 2 * (T0 + soma);
}  // Asymptotics = 4/6 as k -> infinity

template <class Real>
Real serie_impar(int k, Real alpha) {
    Real inv_2k = 2/static_cast<Real>(k);
    Real x = inv_2k*inv_2k;
    Real x_2n = 1;
    Real prod = 1/static_cast<Real>(6);  // (alpha-1)*alpha/(n+1)! para n=2.
    Real T0 = prod;
    Real soma = 0;
    for (int n = 4 ; n < 300; n += 2) {
        prod *= (alpha-3 + n)*(alpha-2 + n)/(n*(n+1));
        x_2n *= x;
        auto soma_old = soma;
        soma = soma_old + prod*x_2n*(3-n);
        if (soma == soma_old) {
            break;
        }
    }
    return (soma + T0)*4;
}  // Asymptotics = 4/6 as k -> infinity

// This function calcultes the coefficient para k = 1 for ealpha != 1.
// The limit for ealpha -> 1 is (8 - 5*std::log(3))/2,
// which comes from:
// (exp(a*x) - 1)/x = a,   x -> 0.
template <class Real>
Real melhorada_1(Real alpha) {  // k = 1 alpha != 1
    Real log3 = std::log(3);
    Real T_1 = -std::expm1(log3*(1 - alpha)) / (1 - alpha);
    Real T_2 = 4 + (alpha + 4) * T_1 / 2;
    return T_2 / (alpha * (2 - alpha));
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

    bool const ealpha_one = b118::almost_equal(ealpha, 1.0);

    {  // {0 < ealpha < 2}

    // Now we treat the general case 0 < ealpha < 2
        bool const ealpha_one = b118::almost_equal(ealpha, 1.0);
        double const ca = ealpha_one ? - frlap::normal_const<1>(1.0)
                                     : - frlap::normal_const<1>(ealpha);

        double const ch = ealpha_one ? + 1.0 / deltax
                                     : + std::pow(deltax, -ealpha);

        // REMARK: using Log-Gamma function and then
        //         exponentiate seems better to mitigate
        //         round-off/approximation errors.
        //         Since lgamma() computes the log of
        //         the abs value of Gamma, we only can
        //         do this because the arguments of
        //         our Gammas are always positive.
        //
        coeffs[0] = ch * M_1_SQRTPI * std::exp2(ealpha) *
            std::exp(std::lgamma(0.5 + 0.5 * ealpha)
                   - std::lgamma(2.0 - 0.5 * ealpha));

        // Obs: coeffs[1] is continuous at ealpha = 1 but we prefered to treat it separately.
    
        if (ealpha_one == true) {  // Could test only ealpha == 1
            coeffs[1] = ca*ch * (8 - 5*std::log(3))/2;
        } else {
            coeffs[1] = ca*ch * melhorada_1(ealpha);
        }

        // The formulation in Taylor series allows us to use
        // the same loop for both cases: alpha = 1 e alpha != 1
        for (std::size_t k = 2; k < n; k+=2)
            coeffs[k] = ca * 2.0 * serie_par(k, ealpha)*std::pow(deltax * k, - ealpha)/k;

        for (std::size_t k = 3; k < n; k+=2)
            coeffs[k] = ca * serie_impar(k, ealpha)*std::pow(deltax * k, - ealpha)/k;
    }
}
