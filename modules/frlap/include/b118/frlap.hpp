/*   libb118
 *
 *   modules/frlap/include/b118/frlap.hpp
 *   
 *   Fractional Laplacian module (main header)
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

#include <cmath>
#include <stdexcept>
#include <b118/numbers.hpp>

namespace b118 {
namespace frlap {

///
/// @brief Function template to compute the normalization constant
///        associated to the singular integral representation of
///        the fractional Laplacian.
///
///        Reference:
///        https://en.wikipedia.org/wiki/Fractional_Laplacian#Singular_Operator
///
/// @tparam Real Floating-point type representing a real number.
/// @tparam n The dimension of the real space `\\mathbb{R}^n`.
///
/// @param ealpha The fractional Laplacian exponent `alpha`. Must lie
///               between `0` and `2` inclusive (be aware that some
///               authors use `s=alpha/2` as the fractional Laplacian
///               exponent instead; these are the two standard notations
///               for the exponent).
///
/// @return The normalization constant.
///
template<unsigned n, typename Real>
inline Real normal_const(Real ealpha) {
    if (std::isless<Real>(ealpha, 0) || std::isgreater<Real>(ealpha, 2))
        throw std::invalid_argument("b118::frlap::normal_const(): "
                                    "ealpha must lie between 0 and 2");
    return ealpha * (static_cast<Real>(2) - ealpha) / static_cast<Real>(2)
                  * std::pow(b118::numbers::inv_sqrtpi_v<Real>,
                             static_cast<Real>(n))
                  * std::exp2(ealpha - 1)
                  * std::exp(std::lgamma((static_cast<Real>(n) + ealpha) /
                                          static_cast<Real>(2))
                           - std::lgamma(static_cast<Real>(2)
                                       - ealpha / static_cast<Real>(2)));
}

}  // end namespace frlap
}  // end namespace b118
