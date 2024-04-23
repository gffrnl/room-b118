/*   libb118
 *
 *   modules/linalg/include/b118/linalg/matrix_kind.hpp
 *   
 *   Kinds of matrices
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

namespace b118 { namespace linalg { namespace matrix_kind {

    template<typename...> struct dense {};
    template<typename...> struct column {};
    template<typename...> struct symmetric_upper {};
    template<typename...> struct symmetric_lower {};
    template<typename...> struct skew_symmetric {};
    template<typename...> struct sparse {};
    template<typename...> struct toeplitz {};
    template<typename...> struct symmetric_toeplitz {};

}}}  // end namespace b118::linalg::matrix_kind
