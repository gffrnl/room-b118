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
#include "./method.hpp"
#include "./generator.hpp"

namespace b118 { namespace frlap { namespace gdm { namespace coefficients {

template<typename Real>
class gorenflo_mainardi : public method<Real>,
                          public generator<Real, gorenflo_mainardi> {
    using generator<Real, gorenflo_mainardi>::ealpha;
    using generator<Real, gorenflo_mainardi>::deltax;

 public:
    gorenflo_mainardi(Real const & ealpha, Real const & deltax)
        : generator<Real, gorenflo_mainardi>(ealpha, deltax),
          is_ealpha_one(b118::almost_equal<Real>(ealpha, 1)),
          ch(is_ealpha_one ? (static_cast<Real>(1)
                                 / (deltax * b118::numbers::pi_v<Real>))
                           : (static_cast<Real>(1)
                                 / pow(deltax, ealpha)
                                 / cos(ealpha * b118::numbers::pi_v<Real>
                                        / static_cast<Real>(2)))),
          c0(is_ealpha_one ? static_cast<Real>(2) : (
            (ealpha < static_cast<Real>(1)) ? static_cast<Real>(1)
                                            : - ealpha)),
          c1((ealpha * (ealpha - static_cast<Real>(1)) + static_cast<Real>(2))
                / static_cast<Real>(4)),
          ck((ealpha < static_cast<Real>(1))     ?
                (ealpha * (ealpha - static_cast<Real>(1))
                        / static_cast<Real>(2))  :
                (ealpha * (static_cast<Real>(2) - ealpha)
                        * (ealpha - static_cast<Real>(1))
                        / static_cast<Real>(2)))
    {}

    Real operator()(std::size_t const & k) const override {
        //   - there is a way to mitigate the large errors
        //     when alpha ~~ 1 but alpha != 1 ???
        //   - when alpha == 1 it is better use 1.0/(k*(k+1))
        //     or take the logs and then exponentiate ???

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
            return ch * c0;

        // REMARK: Using Log-Gamma function and then
        //         exponentiate seems to be better to
        //         mitigate round-off/approximation errors.
        //         We only can do this because the result
        //         of our Gammas are positive, since lgamma()
        //         computes the logarithm of Gamma's ABSOLUTE
        //         value.

        if (is_ealpha_one)
            return - ch / (static_cast<Real>(k) * static_cast<Real>(k + 1));
                // or it is better to use
                // - ch / exp(  log(static_cast<Real>(k))
                //             + log(static_cast<Real>(k+1)) );  ???

        if (ealpha > static_cast<Real>(1)) {
            if (k == 1)
                return ch * c1;

            return ch * ck * exp(lgamma(static_cast<Real>(k + 1) - ealpha)
                               - lgamma(static_cast<Real>(k + 2))
                               - lgamma(static_cast<Real>(3) - ealpha));
        }

        return ch * ck * exp(lgamma(static_cast<Real>(k) - ealpha)
                           - lgamma(static_cast<Real>(k + 1))
                           - lgamma(static_cast<Real>(2) - ealpha));
    }

 private:
    bool const is_ealpha_one;
    Real const ch;
    Real const c0;
    Real const c1;
    Real const ck;
};

}}}}  // end namespace b118::frlap::gdm::coefficients
