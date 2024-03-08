/*   libb118
 *
 *   include/b118/type_traits.hpp
 *   
 *   Type traits
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

namespace b118 {

    // Real type
    // TODO(gffrnl): [ ] Create a struct for this
    //               [ ] Enable expented precision floating-point types
    template<typename T>
    using enable_if_real_t
        = typename std::enable_if<std::is_floating_point<T>::value, T>::type;

}  // end namespace b118
