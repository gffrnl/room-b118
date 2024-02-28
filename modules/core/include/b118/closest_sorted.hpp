/*   libb118
 *
 *   modules/base/b118/almost_equal.hpp
 *
 *   Algorithm that returns the index the element in a sorted range
 *   closest to a given value
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
#include <iterator>

namespace b118 {

    template<class ForwardIt, class T>
    std::size_t closest_sorted(ForwardIt first, ForwardIt last, T value) {
        // Validate arguments
        if (value < *first || value > *std::prev(last))
            throw std::invalid_argument("b118::closest_sorted: "
                                        "value outside the range");

        std::size_t idx =
            std::lower_bound(first, last, value) - first;

        if (idx == 0) return idx;

        ForwardIt const lower = std::next(first, idx - 1);
        ForwardIt const upper = std::next(first, idx);

        // if (idx > 0 && value - *lower < *upper - value) { --idx; }

        if (value - *lower < *upper - value) --idx;

        return idx;
    }

}  // end namespace b118
