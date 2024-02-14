/*   libb118
 *
 *   modules/linalg/b118/fsp.hpp
 * 
 *   Fast symmetric Toeplitz product
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

#include <fftw3.h>
#include <cstddef>
#include <cstring>
#include <cmath>
#include <vector>

std::vector<double> fast_symm_prod(std::vector<double> const& c,
                                   std::vector<double> const& x) {
    // TODO(gffrnl): VALIDATE ARGUMENTS

    std::size_t n = c.size();
    std::size_t k;  // padding
    std::size_t m;  // the dimension of the augmented matrix

    // Augmentation:
    // A1 and x must be embedded into n+2*(k+1) arrays
    std::vector<double> mu;   // first row of augmented matrix
    std::vector<double> y;    // augmented vector
    std::vector<double> aux;  // auxiliary array

    fftw_plan plan;

    // The sizes needed to allocate memory for the augmented arrays
    k = (n % 2 == 0) ? (n - 4) / 2 : (n - 3) / 2;
    m = n + 2 * (k + 1);

    // Augmented arrays allocation
    mu.resize(m, 0);
    y.resize(m, 0);
    aux.resize(m, 0);

    // Construct the augmented arrays
    std::memcpy(mu.data(), c.data(), n * sizeof(double));
    for (size_t i = 2; i < n; ++i)
        mu[i-2] -= c[i];
    std::memcpy(y.data()+k+1, x.data(), n * sizeof(double));

    // Compute the DST1 of mu and store in aux
    plan = fftw_plan_r2r_1d(m, mu.data(), aux.data(),
                            FFTW_RODFT00, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();

    // Compute the DST1 of y and store in mu
    plan = fftw_plan_r2r_1d(m, y.data(), mu.data(),
                            FFTW_RODFT00, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();

    y.resize(0);
    y.shrink_to_fit();  // y no more needed

    // Multiply
    {
        double scaling = 1.0 / (4.0 * (m + 1));
        for (std::size_t i = 0; i < m; ++i)
            aux[i] *= (scaling * mu[i] / std::sin((i+1) * M_PI / (m+1)));
    }
    // NOTE: need to multiply lamb by 1/2 rather
    //         sqrt(n+1)/2 because in fftw3 scaling
    //         of DST1 is 2

    // Now compute the DST1 of xemb
    plan =
    fftw_plan_r2r_1d(m, aux.data(), mu.data(), FFTW_RODFT00, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();

    aux.resize(0);
    aux.shrink_to_fit();  // aux no more needed

    std::vector<double> b(n);

    /* Copy only necessary values to b */
    std::memcpy(b.data(), mu.data()+k+1, n * sizeof(double));

  return b;
}