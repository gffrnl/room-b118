/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/trunc_uniform.hpp
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

#include <cblas.h>
#include <cassert>
#include <type_traits>
#include <b118/frlap/detail/fstp.hpp>
#include <b118/frlap/gdm.hpp>
#include <b118/frlap/gdm/strategy.hpp>
#include <b118/frlap/gdm/strategies/gorenflo_mainardi.hpp>
#include <b118/frlap/gdm/strategies/huang_oberman_linear.hpp>
#include <b118/frlap/gdm/strategies/huang_oberman_quadratic.hpp>
#include <b118/frlap/gdm/strategies/centered_3_point_periodized.hpp>
#include <b118/frlap/gdm/convolucao.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

void conv1d_ingenua(std::size_t output_size, std::size_t input_size,
                    double const * const kernel,
                    double const * const y,
                    double       * const fr);

void cross1d_ingenua(std::size_t output_size, std::size_t input_size,
                     double const * const kernel,
                     double const * const y,
                     double       * const fr);



namespace b118 {
namespace frlap {
namespace gdm {

template<
    class Strategy,
    typename = std::enable_if<
                   std::is_base_of<strategy, Strategy>::value,
                   bool
               >
>
struct trunc_uniform final : public general_differences_method {
    trunc_uniform(double ealpha, double deltax)
        : general_differences_method(ealpha, deltax)
    {}

    trunc_uniform() : trunc_uniform(0.0, 0.0) {}

    Strategy strategy;

    void compute(std::vector<double> const & y     ,
                 std::size_t                ja     ,
                 std::size_t                jb     ,
                 std::vector<double>&       frLap_y) override {
        if (jb > y.size()-1)
        throw "jb > y.size()";
        if (ja > jb)
        throw "ja > jb";

        std::size_t const n  = y.size();
        std::size_t const na = ja;
        std::size_t const nb = n-1-jb;
        std::size_t const n0 = jb-ja+1;
        std::size_t const nc = (jb+1 > n-ja)? jb+1 : n-ja;

        frLap_y.resize(n0);
        frLap_y.shrink_to_fit();

        std::vector<double> mu(nc);
        strategy.generate_coefficients(ealpha, deltax, mu.data(), mu.size());
        if (false) {
            conv1d_ingenua(n0, na, mu.data()+ja-na+1, y.data(), frLap_y.data());

            cross1d_ingenua(n0, nb, mu.data() + 1, y.data() + jb + 1,
                                    frLap_y.data());
        } else {
            std::size_t const conv_size = n0 + std::max(na, nb) - 1;
            convolution conv(conv_size);
            conv.create_plans(conv_size);
            conv.conv(n0, na, mu.data()+ja-na+1, y.data(), frLap_y.data(), false);
            //conv1d_ingenua(n0, na, mu.data()+ja-na+1, y.data(), frLap_y.data());

            conv.conv(n0, nb, mu.data() + 1, y.data() + jb + 1, frLap_y.data(), true);
        }

        {
            std::vector<double> yint(n0);
            // 5.1. Assembly of Yint:
            for (size_t i = 0; i < n0; i++)
                yint[i] = y[i+ja];
            std::vector<double> Ayint(n0);
            // 5.4. Fast symmetric toeplitz-vector A*Yint:
            fast_symm_toeplitz_prod(n0, mu.data(), yint.data(), Ayint.data());
            for (std::size_t i = 0; i < n0; i++)
                frLap_y[i] += Ayint[i];
        }
    }
};

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118


namespace b118 {
namespace frlap {
namespace gdm {

//   typedef trunc_uniform<
//     strategies::spectral
//   > trunc_uniform_spec;

//   typedef trunc_uniform<
//     strategies::spectral_ooura
//   > trunc_uniform_spec_ooura;

//   typedef trunc_uniform<
//     strategies::spectral_qawo
//   > trunc_uniform_spec_qawo;

//   typedef trunc_uniform<
//     strategies::spectral_tanh_sinh
//   > trunc_uniform_spec_thsh;

  using trunc_uniform_gormai =
    trunc_uniform<strategies::gorenflo_mainardi>;

  using trunc_uniform_huob1 =
    trunc_uniform<strategies::huang_oberman_linear>;

  using trunc_uniform_huob2 =
    trunc_uniform<strategies::huang_oberman_quadratic>;

  using trunc_uniform_c3point =
    trunc_uniform<strategies::centered_3_point_periodized>;

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118




// Convolution via dot product
// The output fr satisfies:
// fr[i] += sum(kernel[input_size - 1 - i + j] * y[j], j = 0 ... input_size - 1), i = 0 .. output_size - 1
// kernel[k] must be defined in [0, output_size + input_size - 1]
void conv1d_ingenua(std::size_t output_size, std::size_t input_size,
                    double const * const kernel,
                    double const * const y,
                    double       * const fr) {
    // Obs: Blas shifts the pointer when INCX < 0:
    //  ...
    //     IX = 1
    //     IY = 1
    //     IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    //     IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    //     DO I = 1,N
    //         DTEMP = DTEMP + DX(IX)*DY(IY)
    //         IX = IX + INCX
    //         IY = IY + INCY
    //     END DO
    // ...
    for (std::size_t i = 0; i < output_size; ++i) {
        fr[i] += cblas_ddot(input_size,
                            const_cast<double*>(kernel) + i,
                            -1,      // increment on kernel = -1
                            const_cast<double*>(y),
                            1);
    }
}


// Cross-correlation via dot product
// The output fr satisfies:
// fr[i] += sum(kernel[output_size - 1 - i + j] * y[j], j = 0 ... input_size - 1), i = 0 .. output_size - 1
// kernel[k] must be defined in [0, output_size + input_size - 1]
void cross1d_ingenua(std::size_t output_size, std::size_t input_size,
                     double const * const kernel,
                     double const * const y,
                     double       * const fr) {
    for (std::size_t i = 0; i < output_size; ++i) {
        fr[i] += cblas_ddot(input_size,
                            const_cast<double*>(kernel) + (output_size - 1) - i,
                            1,      // increment on kernel = 1
                            const_cast<double*>(y),
                            1);
    }
}
