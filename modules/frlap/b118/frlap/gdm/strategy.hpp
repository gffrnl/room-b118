/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/strategy.hpp
 *   
 *   Strategy concept to represent a coefficient generation strategy
 *   for a gdm for the fractional Laplacian
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

#include <vector>

namespace b118 {
namespace frlap {
namespace gdm {

struct strategy {
    virtual ~strategy() {}

    virtual
    void generate_coefficients(double  ealpha,
                               double  deltax,
                               double* mu) const = 0;
};

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
