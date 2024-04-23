/*   libb118
 *
 *   modules/base/b118/utility.hpp
 *
 *   Utilities
 *
 *   Copyright (C) 2024  Guilherme F. Fornel        <gffrnl@gmail.com>
 *                       Fabio Souto de Azevedo     <fazedo@gmail.com>
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


template <typename Real>
class equally_spaced {
 public:
    using size_type = std::size_t;

 private:
    Real const a, b;
    size_type const N;
    Real const step;

 public:
    struct iterator {
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type        = Real;
        using difference_type   = std::ptrdiff_t;
        using pointer           = void;
        using reference         = Real &;

        equally_spaced const * const that;
        size_type n;

        iterator& operator++() {
            ++n;
            return *this;
        }

        iterator operator++(int) {
            iterator ret = *this;
            ++(*this);
            return ret;
        }

        iterator& operator--() {
            --n;
            return *this;
        }

        iterator operator--(int) {
            iterator ret = *this;
            --(*this);
            return ret;
        }

        Real operator*() const {
            return (*that)[n];
        }

        bool operator!=(iterator const& rhs) const {
            return n != rhs.n;
        }

        Real operator[] (size_type k) const {
            return (*that)[n + k];
        }
    };

    using const_iterator = iterator const;

 public:
    equally_spaced(Real a, Real b, std::size_t N)
        : a(a), b(b), N(N), step((b-a)/(N-1))
    {}

    size_type size() const {
        return N;
    }

    Real spacing() const {
        return step;
    }

    Real operator[] (size_type n) const {
        Real x = static_cast<Real>(N - 1 - n)/ (N - 1);
        Real y = static_cast<Real>(n)/ (N - 1);
        return x*a + y*b;
        // return a + step*n;
    }

    size_type closest(Real x) {
        auto pos =
            static_cast<size_type>((x - a)/step + static_cast<Real>(1)/2);
        if (pos < 0) {
            return 0;
        } else if (pos >= N) {
            return N;
        } else {
            return pos;
        }
    }

    iterator begin() const {
         return iterator{this, static_cast<size_t>(0)};
    }

    iterator end() const {
        return iterator{this, N};
    }
};

}  // end namespace b118

