/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/spectral_qawo.hpp
 *
 *   Huang & Oberman coefficients - quadratic interpolation
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

#include <gsl/gsl_integration.h>
#include <cmath>
#include <b118/almost_equal.hpp>
#include <b118/numbers.hpp>
#include "../../../frlap.hpp"
#include "./method.hpp"
#include "./generator.hpp"

namespace b118 { namespace frlap { namespace gdm { namespace coefficients {

template<typename Real> struct spectral_qawo;

template<>
class spectral_qawo<double> final : public method<double>,
                                    public generator<double, spectral_qawo> {
    using generator<double, spectral_qawo>::ealpha;
    using generator<double, spectral_qawo>::deltax;

 public:
    spectral_qawo(double const & ealpha, double const & deltax)
        : generator<double, spectral_qawo>(ealpha, deltax),
          c1(1.0 / (b118::numbers::pi_v<double> * deltax)),
          c2(2.0 / (deltax * deltax)),
          ch(1.0 / (b118::numbers::pi_v<double> * std::pow(deltax, ealpha))),
          levels(30) {
        if (!b118::almost_equal<double>(ealpha, 0.0) &&
            !b118::almost_equal<double>(ealpha, 1.0) &&
            !b118::almost_equal<double>(ealpha, 2.0)) {
            levels = 30;
            w = gsl_integration_workspace_alloc(levels);
            t = gsl_integration_qawo_table_alloc(0,
                b118::numbers::pi,
                GSL_INTEG_COSINE, levels);
            expon = this->ealpha;
            F.function = &f;
            F.params = &expon;
            err = 0.0;
        }
    }

    ~spectral_qawo() {
        if (!b118::almost_equal<double>(ealpha, 0.0) &&
            !b118::almost_equal<double>(ealpha, 1.0) &&
            !b118::almost_equal<double>(ealpha, 2.0)) {
            gsl_integration_qawo_table_free(t);
            gsl_integration_workspace_free(w);
        }
    }

    double operator()(std::size_t const & k) const override {
        // The boundary cases ealpha = 0 , ealpha = 1 and ealpha = 2
        if (b118::almost_equal<double>(ealpha, 0.0)) {
            if (k == 0) return 1.0;
            return 0.0;
        }
        if (b118::almost_equal<double>(ealpha, 1.0)) {
            if (k == 0)
                return b118::numbers::pi_v<double> / (2 * deltax);

            if (k%2 == 0)
                return 0.0;

            return - c1 * (2.0 / (k * k));
        }
        if (b118::almost_equal<double>(ealpha, 2.0)) {
            if (k == 0)
                return (b118::numbers::pi_v<double>
                      * b118::numbers::pi_v<double>)
                      / (3.0 * deltax * deltax);

            if (k%2 == 0)
                return c2 / (k * k);

            return - c2 / (k * k);
        }

        // The cases {0 < ealpha < 1} U {1 < ealpha < 2}
        if (k == 0)
            return std::pow(b118::numbers::pi_v<double>/deltax, ealpha)
                / (ealpha + 1.0);

        double value = 0;
        gsl_integration_qawo_table_set(t, k,
            M_PI, GSL_INTEG_COSINE);
        gsl_integration_qawo(&F, 0.0, 1.1e-12, 1.1e-12,
            levels, w, t, &value, &err);
        value *= ch;
        return value;
    }

    double get_err() const { return err; }

 private:
    double const c1;
    double const c2;
    double const ch;

    std::size_t levels;
    mutable double expon;
    mutable gsl_integration_workspace  * w;
    mutable gsl_integration_qawo_table * t;
    mutable gsl_function F;
    mutable double err;
    
    static double f(double x, void* p) {
        double const expon = *static_cast<double *>(p);
        return std::pow(x, expon);
    }
};

}}}}  // end namespace b118::frlap::gdm::coefficients
