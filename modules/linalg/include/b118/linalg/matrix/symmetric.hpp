/*   libb118
 *
 *   modules/linalg/include/b118/linalg/matrix/symmetric.hpp
 *   
 *   Template class to represent symmetric matrices
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
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace b118 { namespace linalg {

    namespace ublas = boost::numeric::ublas;

    template<typename T>
    class matrix<T, matrix_kind::symmetric_lower>
        : public boost::numeric::ublas::symmetric_matrix<
                     T,
                     boost::numeric::ublas::lower
                 > {

     public:  // exports
        using size_type       =
            typename ublas::symmetric_matrix<T, ublas::lower>::size_type;

        using value_type      =
            typename ublas::symmetric_matrix<T, ublas::lower>::value_type;

        using reference       =
            typename ublas::symmetric_matrix<T, ublas::lower>::reference;

        using const_reference =
            typename ublas::symmetric_matrix<T, ublas::lower>::const_reference;

     public:  // constructors
        explicit matrix(size_type sz)
            : ublas::symmetric_matrix<T, ublas::lower>(sz)
        {}

     public:  // accessors
        inline const_reference operator()(size_type i, size_type j) const {
            return
                ublas::symmetric_matrix<T, ublas::lower>::operator()(i-1, j-1);
        }

        inline reference operator()(size_type i, size_type j) {
            return
                ublas::symmetric_matrix<T, ublas::lower>::operator()(i-1, j-1);
        }
    };

}}  // end namespace b118::linalg
