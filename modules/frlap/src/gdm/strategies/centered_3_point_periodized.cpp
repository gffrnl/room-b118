/*   libfrlap
 *
 *   src/gdm/strategies/centered_3_point_periodized.cpp
 *     3 point centered periodized strategy
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

#include <frlap/gdm/strategies/centered_3_point_periodized.hpp>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <cfloat>

static bool almosteq(double a, double b)
{
  double p = std::fabs(a + b);
  double m = std::fabs(a - b);

  if ( p < DBL_EPSILON )
    return true;

  return ( m < (DBL_EPSILON * p) );
}

/**
 *  Math special function beta(a, b);
 *    - here a, b are reals, but in fact can be complex
 *    - definition:
 *        beta(a, b) =
 *          (gamma(a) * gamma(b)) / gamma(a + b) 
 */
static double beta (double a, double b)
{
  // !!! a and b must be positive !!!
  /*
   * REMARK: using log-gamma function and then exponentiate
   *         is better to mitigate round-off/approximation
   *         errors
   */
  return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}


void
frlap::gdm::strategies::centered_3_point_periodized::
generate_coefficients(double                ealpha,
		      double                deltax,
		      std::vector<double> & coeffs) const
{
  const std::size_t n = coeffs.size();

  { /* treat the boundary cases alpha = 0 and alpha = 2 */
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

  { /* {0 < alpha < 2} */
    double ch = 1.0 / ( M_PI * std::pow(deltax, ealpha) );
    double cs = - ch * std::sin(0.5 * ealpha * M_PI);

    coeffs[0] = ch * ( std::exp2(ealpha) * beta(0.5 + 0.5 * ealpha, 0.5) );

    for (std::size_t k = 1; k < n; ++k)
      coeffs[k] = cs * beta((double) k - 0.5 * ealpha, 1.0 + ealpha);
  }
}
