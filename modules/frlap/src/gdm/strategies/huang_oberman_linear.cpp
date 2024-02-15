/*   libfrlap
 *
 *   src/gdm/strategies/huan_oberman_linear.cpp
 *     Huang & Oberman Linear strategy
 *
 *   Copyright (C) 2023  Guilherme F. Fornel <gffrnl@gmail.com>
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

#include <frlap.hpp>
#include <frlap/gdm/strategies/huang_oberman_linear.hpp>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <cfloat>

#ifndef M_1_SQRTPI
#define M_1_SQRTPI 0.564189583547756286948079451561
#endif

static bool almosteq(double a, double b)
{
  double p = std::fabs(a + b);
  double m = std::fabs(a - b);

  if ( p < DBL_EPSILON )
    return true;

  return ( m < (DBL_EPSILON * p) );
}

extern double d1G_alpha_ne_1(double, std::size_t);
extern double d2G_alpha_ne_1(double, std::size_t);
extern double d1G_alpha_eq_1(std::size_t);
extern double d2G_alpha_eq_1(std::size_t);

void
frlap::gdm::strategies::huang_oberman_linear::
generate_coefficients(double                ealpha,
		      double                deltax,
		      std::vector<double> & coeffs) const
{
  const std::size_t n = coeffs.size();

  { // treat the boundary cases alpha = 0 and alpha = 2
    if (almosteq(ealpha, 0.0))
      {
	coeffs[0] = 1.0;
	std::fill(coeffs.begin() + 1, coeffs.end(), 0.0);
	return;
      }

    if (almosteq(ealpha, 2.0))
      {
	coeffs[0] =  2.0 / (deltax * deltax);
	coeffs[1] = -1.0 / (deltax * deltax);
	std::fill(coeffs.begin() + 2, coeffs.end(), 0.0);
	return;
      }
  }

  { // {0 < alpha < 2}

    double ca = 0.0;
    double ch = 0.0;
    
    if (almosteq(ealpha, 1.0))
      {
        ca = - frlap::normal_const<1>(1.0);
        ch = 1.0 / deltax;
      }
    else
      {
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
      exp(  std::lgamma(0.5 * (ealpha + 1.0))
	  - std::lgamma(2.0 -  0.5 *ealpha ) );

    ch *= ca;
    if (almosteq(ealpha, 1.0))
      {
        coeffs[1] = ch * (  1.0
			  - d2G_alpha_eq_1(1)
			  + d1G_alpha_eq_1(2)
			  - d1G_alpha_eq_1(1) );
        for (std::size_t k = 2; k < n; ++k)
          coeffs[k] = ch * (        d1G_alpha_eq_1(k+1)
			    - 2.0 * d1G_alpha_eq_1(k  )
			    +       d1G_alpha_eq_1(k-1) );
      }
    else
      {
        coeffs[1] = ch * (  1.0 / (2.0 - ealpha)
			  - d2G_alpha_ne_1(ealpha, 1)
			  + d1G_alpha_ne_1(ealpha, 2)
			  - d1G_alpha_ne_1(ealpha, 1) );
        for (std::size_t k = 2; k < n; ++k)
          coeffs[k] = ch * (        d1G_alpha_ne_1(ealpha, k+1)
			    - 2.0 * d1G_alpha_ne_1(ealpha, k  )
			    +       d1G_alpha_ne_1(ealpha, k-1) );
      }
  }
}
