/*   libb118
 *
 *   include/b118/linspace.hpp
 *     Contains:
 *       - Linspace algorithms
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
#include <iterator>

namespace b118 {
    enum exclude {
        none  = 0b00,
        left  = 0b01,
        right = 0b10,
        both  = 0b11,
    };
}


namespace b118 {

//
// *** PASSING ITERATORS ***
//

template<
    typename T,
    class ForwardIt
>
T linspace(ForwardIt first, ForwardIt last, T a, T b, char opt = exclude::none) {
    auto n = std::distance(first, last);
    if (n <= 1) throw std::invalid_argument("distance(first, last) <= 1");
    --n;
    a /= static_cast<T>(n);
    b /= static_cast<T>(n);
    if (opt & exclude::both) {
        T const step = (b - a) / static_cast<T>(((opt == exclude::both)? (n + 2) : (n + 1)));
        if (opt & exclude::left ) a += step;
        if (opt & exclude::right) b -= step;
    }
    for (std::size_t k = 0; first != last; ++first, ++k) {
        *first = b * static_cast<T>(k) + a * static_cast<T>(n - k);
    }
    return (b * static_cast<T>(1) + a * static_cast<T>(n - 1)) - (a * static_cast<T>(n));
}

//
// *** PASSING CONTAINERS ***
//

template<
    class T,
    template<class, class...> class Container,
    class... Other
>
T linspace(Container<T, Other...>& x, T a, T b, char opt = exclude::none) {
    return linspace<T>(x.begin(), x.end(), a, b, opt);
}

template<
    class T,
    std::size_t N,
    template<class, std::size_t, class...> class Container,
    class... Other
>
T linspace(Container<T, N, Other...>& x, T a, T b, char opt = exclude::none) {
    return linspace<T>(x.begin(), x.end(), a, b, opt);
}

//
// *** RETURNING CONTAINERS ***
//

template<
    template<typename, class...> class Container,
    typename T,
    class... Other
>
Container<T, Other...> linspace(T a, T b, std::size_t n, char opt = exclude::none, T* step = nullptr) {
    if (n <= 1) throw std::invalid_argument("n <= 1");
    Container<T, Other...> x(n);
    T const aux = linspace<T>(x.begin(), x.end(), a, b, opt);
    if (step != nullptr)
        *step = aux;
    return x;
}

template<
    template<typename, std::size_t, class...> class Container,
    typename T,
    std::size_t N,
    class... Other
>
Container<T, N, Other...> linspace(T a, T b, char opt = exclude::none, T* step = nullptr) {
    static_assert(N > 1);
    Container<T, N, Other...> x;
    T const aux = linspace<T>(x.begin(), x.end(), a, b, opt);
    if (step != nullptr)
        *step = aux;
    return x;
}

} // end namespace b118
