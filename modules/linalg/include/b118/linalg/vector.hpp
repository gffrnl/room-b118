/*   libb118
 *
 *   modules/linalg/include/b118/linalg/vector.hpp
 *   
 *   b118::linalg::vector is a syntatic sugar for a column matrix
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

#include "./matrix_kind.hpp"
#include "./matrix/column.hpp"

namespace b118 { namespace linalg {

    template<typename T>
    using vector = matrix<T, matrix_kind::column>;

}}  // end namespace b118::linalg
