/*   libb118
 *
 *   modules/linalg/include/b118/linalg/toeplitz/solvers/durbin.hpp
 *   
 *   TODO(gffrnl): DESCRIPTION
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

#include <iterator>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace b118 { namespace linalg { namespace toeplitz { namespace solvers {

template<
    typename Real,
    class InputBidirIt,
    class OutputBidirIt
>
OutputBidirIt durbin(InputBidirIt  rhs_beg, InputBidirIt rhs_end,
                     OutputBidirIt sol_beg) {
    {
        auto const len = std::distance(rhs_beg, rhs_end);

        if (len < static_cast<decltype(len)>(0))
            throw std::invalid_argument("durbin(): "
                                        "std::distance(rhs_beg, rhs_end) < 0");

        if (len == static_cast<decltype(len)>(0))
            return sol_beg;

        if (len == static_cast<decltype(len)>(1)) {
            *sol_beg = *rhs_beg;
            return std::next(sol_beg);
        }
    }


    auto rhs_it = rhs_beg;
    auto sol_it = std::next(sol_beg);
    auto sol_rend = std::make_reverse_iterator(sol_beg);
    auto sol_rit  = sol_rend;

    Real alpha, beta;

    alpha = *sol_beg = - (*rhs_it);
    beta = static_cast<Real>(1);

    auto op = [&alpha](Real const & a, Real const & b) -> Real {
                          return alpha * a + b;
                      };

    while (++rhs_it != rhs_end) {
        ++sol_it;
        --sol_rit;

        beta  *= (static_cast<Real>(1) - alpha * alpha);

        alpha = - (*rhs_it + std::transform_reduce(rhs_beg, rhs_it, sol_rit,
                                                   static_cast<Real>(0)))
            / beta;

        std::transform(sol_rit, sol_rend, sol_beg, std::prev(sol_rit),
                       op);
        *sol_beg = alpha;
        std::reverse(sol_beg, sol_it);
    }

    return ++sol_it;
}

}}}}  // end namespace b118::linalg::toeplitz::solvers
