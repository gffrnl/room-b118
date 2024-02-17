/*   libb118
 *
 *   modules/base/b118/utility.hpp
 *
 *   Utilities
 *
 *   Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>
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

namespace b118 {

// TODO(gffrnl,fazedo) DESCRIPTION
// TODO(gffrnl,fazedo) BETTER IMPLEMENTATION ??
template<typename Unsigned>
Unsigned next_exp2(Unsigned n) {
    Unsigned x = 1;
    while (x < n) {
        x <<= 1;
    }
    return x;
}

}  // end namespace b118
