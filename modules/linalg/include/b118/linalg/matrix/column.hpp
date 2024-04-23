/*   libb118
 *
 *   modules/linalg/include/b118/linalg/matrix/column.hpp
 *   
 *   Template class to represent column matrices
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

#include "./matrix.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace b118 { namespace linalg {

    namespace ublas = boost::numeric::ublas;

    template<typename T>
    class matrix<T, matrix_kind::column>
        : public boost::numeric::ublas::vector<T> {

     public:  // exports
        using size_type       = typename ublas::vector<T>::size_type;
        using value_type      = typename ublas::vector<T>::value_type;
        using reference       = typename ublas::vector<T>::reference;
        using const_reference = typename ublas::vector<T>::const_reference;

     public:  // constructors
        explicit matrix(size_type sz) : ublas::vector<T>(sz) {}

     public:
        using ublas::vector<T>::size;

     public:  // accessors
        inline const_reference operator()(size_type i) const {
            return ublas::vector<T>::operator()(i-1);
        }

        inline reference operator()(size_type i) {
            return ublas::vector<T>::operator()(i-1);
        }

        inline const_reference operator[](size_type i) const {
            return ublas::vector<T>::operator()(i);
        }

        inline reference operator[](size_type i) {
            return ublas::vector<T>::operator()(i);
        }
    };

}}  // end namespace b118::linalg