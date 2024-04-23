/*   libb118
 *
 *   modules/linalg/include/b118/linalg/toeplitz/solvers/trench.hpp
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
#include "./durbin.hpp"

namespace b118 { namespace linalg { namespace toeplitz { namespace solvers {

template<
    typename Real,
    class InputBidirIt,
    class OutputBidirIt
>
OutputBidirIt trench(InputBidirIt  row_beg, InputBidirIt row_end,
                     OutputBidirIt inv_beg) {
    std::size_t n = 0;  // system's size

    {
        auto const len = std::distance(row_beg, row_end);

        if (len < static_cast<decltype(len)>(0))
            throw std::invalid_argument("trench(): "
                                        "std::distance(row_beg, row_end) < 0");

        if (len == static_cast<decltype(len)>(0))
            return inv_beg;

        if (std::abs(*row_beg) < std::numeric_limits<Real>::epsilon())
            throw std::invalid_argument("trench(): zero diagonal");

        if (len == static_cast<decltype(len)>(1)) {
            *inv_beg = static_cast<Real>(1) / *row_beg;
            return std::next(inv_beg);
        }

        n = static_cast<std::size_t>(len);
    }

    Real const diag = *row_beg;

    if (diag != 1) {  // scale if not one in diagonal
        std::for_each(row_beg, row_end, [&diag](Real & x) { x /= diag; });
    }

    std::vector<Real> y(n-1);  // auxiliary buffer

    durbin<double>(std::next(row_beg), row_end, y.begin());

    Real const gamma = static_cast<Real>(1) / (
        static_cast<Real>(1) + std::transform_reduce(y.cbegin(), y.cend(),
                                                     std::next(row_beg),
                                                     static_cast<Real>(0))
    );


    auto const inv_beg0 = inv_beg;

    *inv_beg = gamma;

    std::transform(y.cbegin(), y.cend(), std::next(inv_beg),
                   [&gamma](Real const & x) -> Real { return gamma * x; });

    std::ptrdiff_t const imax = (n - 1) / 2 + 1;  // it floors the number
    for (std::ptrdiff_t i = 2; i <= imax; ++i) {
        auto inv_it = std::next(inv_beg, n - 2 * (i - 2));
        for (std::ptrdiff_t j = i; j <= n - i + 1; ++j)
            *(inv_it + (j - i)) = *(inv_beg++) + gamma * (
                  y.at(i - 2) * y.at(j - 2) - y.at(n - j) * y.at(n - i) 
            );
        inv_beg = inv_it;
    }

    if (diag != 1) {  // scale if not one in diagonal
        std::for_each(inv_beg0, std::next(inv_beg0, n),
                      [&diag](Real & x) { x /= diag; });
    }

    return (n % 2 == 0) ? (inv_beg += 2) : ++inv_beg;
}

}}}}  // end namespace b118::linalg::toeplitz::solvers
