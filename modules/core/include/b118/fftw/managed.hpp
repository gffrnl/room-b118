/*   libb118
 *
 *   modules/base/b118/fftw/managed.hpp
 *
 *   Base class for fftw3 managed objects
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

#include "./manager.hpp"

#include <iostream>

namespace b118 {
namespace fftw {

// template<typename Real>
// class managed {
//     unsigned num_users;

//  protected:
//     managed() {
//         static unsigned instances = 0;
//         num_users = ++instances;
//     }

//     virtual ~managed() {
//         if (--num_users == 0) b118::fftw::manager<Real>::cleanup();
//         std::cout << "fftw_managed users = " << num_users << std::endl;
//     }

//  public:
//     unsigned users() const { return num_users; }
// };

template<typename Real>
class managed {
    static unsigned instances;

 protected:
    managed() {
        ++instances;
    }

    virtual ~managed() {
        if (--instances == 0) b118::fftw::manager<Real>::cleanup();
    }

 public:
    static unsigned num_users() { return instances; }
};

template<typename Real>
unsigned managed<Real>::instances = 0;

}  // end namespace fftw
}  // end namespace b118
