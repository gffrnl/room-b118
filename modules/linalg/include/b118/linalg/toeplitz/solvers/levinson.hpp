/*   libb118
 *
 *   modules/linalg/include/b118/linalg/toeplitz/solvers/levinson.hpp
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

#include <cstddef>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <limits>

namespace b118 { namespace linalg { namespace toeplitz { namespace solvers {

template<
    typename Real,
    class InputBidirIt1, class InputBidirIt2,
    class OutputBidirIt
>
OutputBidirIt levinson(InputBidirIt1 row_beg, InputBidirIt1 row_end,
                       InputBidirIt2 rhs_beg,
                       OutputBidirIt sol_beg) {
    std::size_t n = 0;  // system's size

    {
        auto const len = std::distance(row_beg, row_end);

        if (len < static_cast<decltype(len)>(0))
            throw std::invalid_argument("levinson(): "
                                        "std::distance(row_beg, row_end) < 0");

        if (len == static_cast<decltype(len)>(0))
            return sol_beg;

        if (std::abs(*row_beg) < std::numeric_limits<Real>::epsilon())
            throw std::invalid_argument("levinson(): zero diagonal");

        if (len == static_cast<decltype(len)>(1)) {
            *sol_beg = *rhs_beg / *row_beg;
            return std::next(sol_beg);
        }

        n = static_cast<std::size_t>(len);
    }

    std::vector<Real> aux(n);  // auxiliary buffer

    auto row_it   = std::next(row_beg);
    auto row_rit  = --std::make_reverse_iterator(row_beg);
    auto rhs_it   = rhs_beg;
    auto sol_it   = sol_beg;
    auto aux_beg  = aux.begin();
    auto aux_it   = aux_beg;
    auto aux_rend = std::make_reverse_iterator(aux_beg);
    auto aux_rit  = aux_rend;

    Real alpha, beta, mu;

    alpha = *aux_beg = - (*row_it / *row_beg);
    beta  = static_cast<Real>(1);
    *sol_beg = (*rhs_beg / *row_beg);

    auto op_mu =
           [&mu](Real const & a, Real const & b) -> Real {
                    return a + mu * b;
                };
    auto op_alpha =
        [&alpha](Real const & a, Real const & b) -> Real {
                    return alpha * a + b;
                };

    while (row_it != row_end) {
        ++row_it;
        --row_rit;
        ++rhs_it;
        ++sol_it;
        ++aux_it;
        --aux_rit;

        beta *= (static_cast<Real>(1) - alpha * alpha);

        mu = (*rhs_it - std::transform_reduce(sol_beg, sol_it, row_rit,
                                              static_cast<Real>(0)))
            / beta / *row_beg;

        std::transform(sol_beg, sol_it, aux_rit, sol_beg, op_mu);
        *sol_it = mu;

        if (row_it != row_end) {
            alpha = - (*row_it + std::transform_reduce(aux_beg, aux_it, row_rit,
                                                       static_cast<Real>(0)))
                / beta / *row_beg;

            std::transform(aux_rit, aux_rend, aux_beg, std::prev(aux_rit),
                           op_alpha);
            *aux_beg = alpha;
            std::reverse(aux_beg, std::next(aux_it));
        }
    }

    return ++sol_it;
}

}}}}  // end namespace b118::linalg::toeplitz::solvers
