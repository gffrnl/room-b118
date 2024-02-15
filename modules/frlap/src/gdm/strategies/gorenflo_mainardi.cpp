/*   libfrlap
 *
 *   src/gdm/strategies/gorenflo_mainardi.cpp
 *     Gorenflo & Mainardi strategy
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

#include <frlap/gdm/strategies/gorenflo_mainardi.hpp>
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

void
frlap::gdm::strategies::gorenflo_mainardi::
generate_coefficients(double                ealpha,
		      double                deltax,
		      std::vector<double> & coeffs) const
{
  // TODO:
  //   - there is a way to mitigate the large errors
  //     when alpha ~~ 1 but alpha != 1 ???
  //   - when alpha == 1 it is better use 1.0/(k*(k+1))
  //     or take the logs and then exponentiate ???

  const std::size_t n = coeffs.size();

  { // treat the boundary cases ealpha = 0 and ealpha = 2
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

  { // {0 < ealpha < 2}
    if (almosteq(ealpha, 1.0))
      {
        double c1 = 1.0 / (deltax * M_PI);
        
        coeffs[0] = c1 * 2.0;

        for (std::size_t k = 1; k < n; ++k)
          coeffs[k] = - c1 / (k * (k+1));
	/* or it is better to use ??
	   coeffs[k] = - c1 / std::exp(  std::log((double)  k   )
	   + std::log((double) (k+1)) );
	*/
      }
    else
      {
        double ch = 1.0 / std::pow(deltax, ealpha);
        double cc = 1.0 / cos(0.5 * ealpha * M_PI);
        double ca = ealpha;

        // REMARK: Using Log-Gamma function and then
        //         exponentiate seems to be better to
        //         mitigate round-off/approximation errors.
        //         We only can do this because the result
        //         of our Gammas are positive, since lgamma()
        //         computes the logarithm of Gamma's ABSOLUTE
        //         value.

        if (ealpha < 1.0)
          {
            coeffs[0] = cc * ch;
            
            ca *= 0.5 * (ealpha - 1.0) * ch;
            for (std::size_t k = 1; k < n; k++)
              coeffs[k] = cc * ( ca *
	          std::exp(  std::lgamma((double) k - ealpha)
			   - std::lgamma(k + 1)
			   - std::lgamma(2.0 - ealpha)        ) );
          }
        else // 1 < ealpha < 2
          {
            coeffs[0] = - cc * (ca * ch);

            ca *= 0.5 * (ealpha - 1.0);
            coeffs[1] = cc * ( (0.5 + 0.5 * ca) * ch);

            ca *= (2.0 - ealpha) * ch;
            for (std::size_t k = 2; k < n; k++)
              coeffs[k] = cc * ( ca *
	          std::exp(  std::lgamma((double) (k+1) - ealpha)
			   - std::lgamma(k + 2)
			   - std::lgamma(3.0 - ealpha)          ) );
          }
      }
  }
}
