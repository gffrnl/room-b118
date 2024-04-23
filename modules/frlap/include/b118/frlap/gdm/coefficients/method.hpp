/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/method.hpp
 *
 *   Some description
 *
 *   Copyright (C) 2024   Guilherme F. Fornel        <gffrnl@gmail.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as publeished by
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

namespace b118 { namespace frlap { namespace gdm { namespace coefficients {

template<typename Real>
struct method {
    virtual ~method() {}
    virtual Real operator()(std::size_t const & k) const = 0;
};

}}}}  // end namespace b118::frlap::gdm::coefficients
