/*   libb118
 *
 *   modules/frlap/src/gdm/strategies/huang_oberman.cpp
 *   
 *   Common functions used by Huang & Oberman strategies
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

#include <cstddef>
#include <cmath>

double d0G_alpha_ne_1(double ealpha, std::size_t k) {
    return std::pow(static_cast<double>(k), 2.0 - ealpha)
        / (ealpha * (ealpha - 1.0) * (2.0 - ealpha));
}

double d1G_alpha_ne_1(double ealpha, std::size_t k) {
    return std::pow(static_cast<double>(k), 1.0 - ealpha)
        / (ealpha * (ealpha - 1.0));
}

double d2G_alpha_ne_1(double ealpha, std::size_t k) {
    return - std::pow(static_cast<double>(k), -ealpha) / ealpha;
}

double d0G_alpha_eq_1(std::size_t k) {
    return static_cast<double>(k) - std::log(static_cast<double>(k)) * k;
}

double d1G_alpha_eq_1(std::size_t k) {
    return - std::log(static_cast<double>(k));
}

double d2G_alpha_eq_1(std::size_t k) {
    return - 1.0 / static_cast<double>(k);
}
