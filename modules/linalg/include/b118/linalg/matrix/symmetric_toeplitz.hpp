/*   libb118
 *
 *   modules/linalg/include/b118/linalg/matrix/symmetric_toeplitz.hpp
 *   
 *   Template class to represent symmetric Toeplitz matrices
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

#include <vector>
#include <memory>
#include <stdexcept>
#include "./matrix.hpp"
#include "./column.hpp"
#include "../toeplitz/fast_symm_toeplitz_prod.hpp"

namespace b118 { namespace linalg {

    template<typename T>
    class matrix<T, matrix_kind::symmetric_toeplitz> {
        std::vector<T> data;

     public:  // exports
        // using size_type       = typename std::size_t;
        // using value_type      = T;
        // using reference       = T&;
        // using const_reference = T const&;

        using size_type       = typename std::vector<T>::size_type;
        using value_type      = typename std::vector<T>::value_type;
        using reference       = typename std::vector<T>::reference;
        using const_reference = typename std::vector<T>::const_reference;

        using iterator        = typename std::vector<T>::iterator;
        using const_iterator  = typename std::vector<T>::const_iterator;

     public:  // constructors
        explicit matrix(size_type sz) : data(sz) {}
        explicit matrix(std::initializer_list<T> init) : data(init) {}

     public:  // destructor
        virtual ~matrix() {}

     public:  // accessors
        inline const_reference operator()(size_type i, size_type j) const {
    #ifdef NDEBUG
            return data[(i > j)? i-j : j-i];
    #else
            return data.at((i > j)? i-j : j-i);
    #endif
        }
        inline reference operator()(size_type i, size_type j) {
    #ifdef NDEBUG
            return data[(i > j)? i-j : j-i];
    #else
            return data.at((i > j)? i-j : j-i);
    #endif
        }

     public:
        inline size_type size1() const { return data.size(); }
        inline size_type size2() const { return size1(); }

        inline iterator begin() { return data.begin(); }
        inline iterator end() { return data.end(); }
        inline const_iterator cbegin() const { return data.cbegin(); }
        inline const_iterator cend() const { return data.cend(); }

    //  public:  // TODO(gffrnl): REMOVE
    //     std::vector<T> get_data() { return data; }
    };

}}  // end namespace b118::linalg


template<typename Real>
b118::linalg::matrix<Real, b118::linalg::matrix_kind::column>
    operator*(b118::linalg::matrix<
                  Real,
                  b118::linalg::matrix_kind::symmetric_toeplitz
              > const & A,
              b118::linalg::matrix<
                  Real,
                  b118::linalg::matrix_kind::column
              > const & x) {
        namespace linalg   = b118::linalg;
        namespace toeplitz = b118::linalg::toeplitz;

        if (x.size() != A.size2())
            throw std::invalid_argument("A and x are not conform in sizes");

        linalg::matrix<Real, linalg::matrix_kind::column> y(x.size());

        toeplitz::fast_symm_toeplitz_prod<Real>{}(A.cbegin(), A.cend(),
                                                  x.cbegin(), y.begin());
        return y;
}


template<class E, class T, class U>
std::basic_ostream<E, T> & operator<<(
    std::basic_ostream<E, T> & os,
    b118::linalg::matrix<U,
                         b118::linalg::matrix_kind::symmetric_toeplitz
    > const & m) {
    using size_type = typename
        b118::linalg::matrix<
            U,
            b118::linalg::matrix_kind::symmetric_toeplitz
        >::size_type;

    size_type size1 = m.size1();
    size_type size2 = m.size2();

    std::basic_ostringstream<E, T, std::allocator<E> > s;

    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision());
    s << '[' << size1 << ',' << size2 << "](";

    if (size1 > 0) {
        s << '(';
        if (size2 > 0)
            s << m(0, 0);
        for (size_type j = 1; j < size2; ++j)
            s << ',' << m(0, j);
        s << ')';
    }

    for (size_type i = 1; i < size1; ++i) {
        s << ",(";
        if (size2 > 0)
            s << m(i, 0);
        for (size_type j = 1; j < size2; ++j)
            s << ',' << m(i, j);
        s << ')';
    }
    s << ')';

    return os << s.str().c_str();
}
