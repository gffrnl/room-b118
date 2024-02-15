/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm.hpp
 *   
 *   Fractional Laplacian module (main header)
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

struct general_differences_method {
    double ealpha;
    double deltax;
    bool verbose;

    general_differences_method(double frac_expon, double grid_step)
        : ealpha{frac_expon}, deltax{grid_step}, verbose{false}
    {}

    virtual ~general_differences_method() {}

    virtual void compute(std::vector<double> const & y     ,
                         std::size_t                ja     ,
                         std::size_t                jb     ,
                         std::vector<double>&       frLap_y) = 0;
};

// Syntatic sugar
using gdm_t = struct general_differences_method;

}  // end namespace frlap
}  // end namespace b118
