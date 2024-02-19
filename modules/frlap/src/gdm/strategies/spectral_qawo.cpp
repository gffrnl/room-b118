/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/spectral_qawo.cpp
 *
 *   Spectral strategy (uses aqwo integrator)
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

#include <gsl/gsl_integration.h>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <b118/almost_equal.hpp>
#include <b118/frlap/gdm/strategies/spectral_qawo.hpp>

static double f(double x, void* p) {
  double const expon = *static_cast<double *>(p);
  return std::pow(x, expon);
}

using b118::frlap::gdm::strategies::spectral_qawo;

void spectral_qawo::generate_coefficients(double      ealpha,
                                          double      deltax,
                                          double*     coeffs,
                                          std::size_t n) const {
  {  // Treat the boundary cases ealpha = 0 and ealpha = 2
    if (b118::almost_equal(ealpha, 0.0)) {
      coeffs[0] = 1.0;
      std::fill(coeffs + 1, coeffs + n, 0.0);
      return;
    }

    if (b118::almost_equal(ealpha, 2.0)) {
      double const c2 = 2.0 / (deltax * deltax);
      coeffs[0] =  (M_PI * M_PI) / (3.0 * deltax * deltax);
      for (std::size_t k = 1; k < n; k+=2)
        coeffs[k] = - c2 / (k * k);
      for (std::size_t k = 2; k < n; k+=2)
        coeffs[k] =   c2 / (k * k);
      return;
    }
  }

  { // Treat the case ealpha = 1
    if (b118::almost_equal(ealpha, 1.0)) {
      double const c1 = 1.0 / (M_PI * deltax);
      coeffs[0] =  M_PI / (2 * deltax);
      for (std::size_t k = 1; k < n; k+=2)
        coeffs[k] = - c1 * (2.0 / (k * k));
      for (std::size_t k = 2; k < n; k+=2)
        coeffs[k] = 0.0;
      return;
    }
  }

  { // Treat the cases {0 < ealpha < 1} U {1 < ealpha < 2}
    double err;
    double ch = 1.0 / (M_PI * std::pow(deltax, ealpha));

    coeffs[0] = std::pow(M_PI/deltax, ealpha) / (ealpha + 1.0);

    {
      gsl_integration_workspace *w;
      gsl_integration_qawo_table* t;
      gsl_function F;
      size_t levels = 30;

      F.function = &f;
      F.params = &ealpha;

      w = gsl_integration_workspace_alloc(levels);
      t = gsl_integration_qawo_table_alloc(0, M_PI, GSL_INTEG_COSINE, levels);

      for (std::size_t k = 1; k < n; ++k) {
        gsl_integration_qawo_table_set(t, k, M_PI, GSL_INTEG_COSINE);
        gsl_integration_qawo(&F, 0.0, 1.1e-12, 1.1e-12, levels, w, t,
                             &coeffs[k],
                             &err);
        coeffs[k] *= ch;
      }

      gsl_integration_qawo_table_free(t);
      gsl_integration_workspace_free(w);
    }
  }
}
