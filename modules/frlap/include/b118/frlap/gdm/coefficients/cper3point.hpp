/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/cper3point.hpp
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

#include <b118/almost_equal.hpp>
#include <b118/numbers.hpp>
#include "./method.hpp"
#include "./generator.hpp"

namespace b118 { namespace frlap { namespace gdm { namespace coefficients {

template<typename Real>
class cper3point final : public method<Real>,
                         public generator<Real, cper3point> {
    using generator<Real, cper3point>::ealpha;
    using generator<Real, cper3point>::deltax;

 public:
    cper3point(Real const & ealpha, Real const & deltax)
        : generator<Real, cper3point>(ealpha, deltax),
          ch(static_cast<Real>(1) /
                (b118::numbers::pi_v<Real> * pow(deltax, ealpha))),
          cs(- ch * sin(ealpha * b118::numbers::pi_v<Real>
                / static_cast<Real>(2)))
    {}

    Real operator()(std::size_t const & k) const override {
        // The boundary cases ealpha = 0 and ealpha = 2
        if (b118::almost_equal<Real>(ealpha, 0)) {
            if (k == 0) return static_cast<Real>(1);
            return static_cast<Real>(0);
        }
        if (b118::almost_equal<Real>(ealpha, 2)) {
            if (k == 0)
                return   static_cast<Real>(2) / (deltax * deltax);
            if (k == 1)
                return - static_cast<Real>(1) / (deltax * deltax);
            return static_cast<Real>(0);
        }

        // The general case 0 < ealpha < 2
        if (k == 0)
            return ch * (exp2(ealpha)
                      * std::beta((ealpha + 1) / static_cast<Real>(2),
                                  static_cast<Real>(1) / static_cast<Real>(2)));
        return cs * std::beta(static_cast<Real>(k)
                                - ealpha / static_cast<Real>(2),
                              ealpha + static_cast<Real>(1));
    }

 private:
    Real const ch;
    Real const cs;
};

}}}}  // end namespace b118::frlap::gdm::coefficients
