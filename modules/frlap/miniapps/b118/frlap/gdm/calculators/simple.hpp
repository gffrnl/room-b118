/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/calculators/simple.hpp
 *
 *   Simple fractional Laplacian calculator using general differences
 *   approximations
 *
 *   Copyright (C) 2023  Guilherme F. Fornel <gffrnl@gmail.com>
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

#include <string>
#include <vector>
#include <memory>
#include <b118/frlap/gdm.hpp>

namespace b118  {
namespace frlap {
namespace gdm   {
namespace calculators  {

struct Interval {
    double a;
    double b;
};

struct Problem {
    virtual ~Problem() {}

    std::string label = "unknown";
    double      ealpha;  // fractional Laplacian exponent alpha
    Interval    domain;  // computational domain
    Interval    viswin;  // visualization domain (viswindow $\subset$ domain)
    std::size_t numnod;  // number of nodes

    std::unique_ptr<b118::frlap::gdm_t> method_ptr;
    // the general differences method

    virtual double y(double x) const = 0;
    virtual double frLap_y(double x) const = 0;
    virtual void compute_far_field(std::vector<double> const& x ,
                                   std::size_t                ja,
                                   std::size_t                jb,
                                   double*                    ff) const = 0;
};

extern std::vector<std::unique_ptr<Problem>> prob_ptrs;

}  // end namespace calculators
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
