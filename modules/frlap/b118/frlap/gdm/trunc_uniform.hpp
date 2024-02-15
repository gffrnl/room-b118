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
#include <b118/frlap/detail/toep.h>
#include <cassert>
#include <type_traits>
#include <b118/frlap/gdm.hpp>
#include <b118/frlap/gdm/strategy.hpp>
#include <b118/frlap/gdm/strategies/huang_oberman_quadratic.hpp>


#ifdef FRLAP_MULTIPLICA_FABIO_BLAS
#include <iostream>
void multiplica_a(size_t n0, size_t na, int ja,
                const std::vector<double>& mu,
                const std::vector<double>& y,
                std::vector<double> *fr);
#endif
#ifdef FRLAP_CONVCORR_INGENUA
#include <iostream>
#include <vector>
// void conv1d_ingenua(size_t n0, size_t na,
//                                     int ja,
//                                     std::vector<double> const & mu,
//                                     std::vector<double> const & y,
//                                     std::vector<double> *fr);

void conv1d_ingenua(size_t n0, size_t na,
                    double const * const mu,
                    double const * const y,
                    double       * const fr);

void cross1d_ingenua(size_t n0, size_t nb,
                                       int jb,
                                       std::vector<double> const & mu,
                                       std::vector<double> const & y,
                                       std::vector<double> *fr);



#endif

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

        {
#ifdef FRLAP_MULTIPLICA_FABIO_BLAS
        std::cout << "trunc_uniform(): "
                  << "using multiplica_a()" << std::endl;
        multiplica_a(n0, na, ja, mu, y, &frLap_y);
#elif defined(FRLAP_CONVCORR_INGENUA)
        std::cout << "trunc_uniform(): "
                  << "conv1d_ingenua()" << std::endl;
        conv1d_ingenua(n0, na, mu.data()+ja, y.data(), frLap_y.data());
#else
        std::vector<double> Ba(n0 * na);
        std::vector<double> Ya(na);
        for (std::size_t i = 0; i < n0; ++i)
            for (std::ptrdiff_t j = 0; j < na; ++j)
            Ba[i*na + j]=mu[ja+i-j];
        // y = ealpha*A*x + beta*y  com ealpha = 1.0 e beta = 0.0
        cblas_dgemv(CblasRowMajor,  // Armazenamento por linha
            CblasNoTrans,   // Não transpor A
            n0, na,         // Dimensões de A
            1.0,            // ealpha
            Ba.data(), // Matriz A
            na,             // LDA specifies the first dimension of A as declared in the calling (sub) program.
            y.data(),             // Vetor x
            1,              // Passo de x
            0.0,            // beta
            frLap_y.data(), // Vetor y
            1);             // Passo de y
#endif
        }

        {
#if defined(FRLAP_CONVCORR_INGENUA)
            {
                std::vector<double> cross_mu_y(n0);
                std::cout << "trunc_uniform(): "
                        << "cross1d_ingenua()" << std::endl;
                cross1d_ingenua(n0, nb, jb, mu, y, &cross_mu_y);
                for (std::size_t i = 0; i < n0; ++i)
                    frLap_y.at(i) += cross_mu_y.at(i);
            }
#else
        std::vector<double> Bb(n0 * nb);
        for (std::ptrdiff_t i = 0; i < n0; ++i)
            for (std::size_t j = 0; j < nb; ++j)
            Bb[i*nb + j]=mu[n0-i+j];
        // y = ealpha*A*x + beta*y  com ealpha = 1.0 e beta = 1.0
        cblas_dgemv(CblasRowMajor,  // Armazenamento por linha
            CblasNoTrans,   // Não transpor A
            n0, nb,         // Dimensões de A
            1.0,            // ealpha
            Bb.data(), // Matriz A
            na,             // LDA specifies the first dimension of A as declared in the calling (sub) program.
            y.data()+jb+1,  // Vetor x
            1,              // Passo de x
            1.0,            // beta
            frLap_y.data(), // Vetor y
            1);             // Passo de y
#endif
        }

        {
        std::vector<double> yint(n0);
        // 5.1. Assembly of Yint:
        for (size_t i = 0; i < n0; i++)
            yint[i] = y[i+ja];
        std::vector<double> Ayint(n0);
        // 5.4. Fast symmetric toeplitz-vector A*Yint:
        int ret = fast_symm_toeplitz_prod(n0, mu.data(), yint.data(), Ayint.data());
        assert(ret == 0); // TODO: change to an error
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
    
//   typedef trunc_uniform<
//     strategies::gorenflo_mainardi
//   > trunc_uniform_gormai;

//   typedef trunc_uniform<
//     strategies::huang_oberman_linear
//   > trunc_uniform_huob1;

//   typedef trunc_uniform<
//     strategies::huang_oberman_quadratic
//   > trunc_uniform_huob2;

    using trunc_uniform_huob2 =
        trunc_uniform<strategies::huang_oberman_quadratic>;

//   typedef trunc_uniform<
//     strategies::centered_3_point_periodized
//   > trunc_uniform_3point;

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118


#ifdef FRLAP_MULTIPLICA_FABIO_BLAS
// conv1d_naïf
void multiplica_a(size_t n0, size_t na, int ja,
                const std::vector<double>& mu,
                const std::vector<double>& y,
                std::vector<double> *fr) {

    //for (std::size_t i = 0; i < n0; ++i)
    //    for (std::ptrdiff_t j = 0; j < na; ++j)
    //        Ba[i*na + j] = mu[ja+i-j];

    // frLap_y[i] = produto interno da linha i de Ba com y
    // a linha i de Ba é dada por mu[ja+i-j] com j variando de 0 a na-1

    // std::vector<double> frLap_y(n0);
    for (std::size_t i = 0; i < n0; ++i) {
        (*fr)[i] = cblas_ddot(
                        na,
                        mu.data() + ja + i + (-na + 1),
                        -1,      // incremento em mu = -1
                        y.data(),
                        1);
    }
    // Aqui tem um gotcha: se o incremento é negativo, a soma começa do final:
    /*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
    */
}
#endif

#ifdef FRLAP_CONVCORR_INGENUA
    void conv1d_ingenua(size_t n0, size_t na,
                        double const * const mu,
                        double const * const y,
                        double       * const fr) {

        // for (std::size_t i = 0; i < n0; ++i)
        //     for (std::ptrdiff_t j = 0; j < na; ++j)
        //         Ba[i*na + j] = mu[ja+i-j];

        // conv_y[i] = produto interno da linha i de Ba com y
        // a linha i de Ba é dada por mu[ja+i-j] com j variando de 0 a na-1

        // std::vector<double> conv_y(n0);
        for (std::size_t i = 0; i < n0; ++i) {
            fr[i] += cblas_ddot(
                            na,
                            //const_cast<double*>(mu) + ja + i + (-na + 1),
                            const_cast<double*>(mu) + i + (-na + 1),
                            -1,      // incremento em mu = -1
                            const_cast<double*>(y),
                            1);
        }
        // Aqui tem um gotcha: se o incremento é negativo, a soma começa
        // do final:
        //
        //     IX = 1
        //     IY = 1
        //     IF (INCX.LT.0) IX = (-N+1)*INCX + 1
        //     IF (INCY.LT.0) IY = (-N+1)*INCY + 1
        //     DO I = 1,N
        //         DTEMP = DTEMP + DX(IX)*DY(IY)
        //         IX = IX + INCX
        //         IY = IY + INCY
        //     END DO
        // END IF
    }

    void cross1d_ingenua(size_t n0, size_t nb,
                                       int jb,
                                       std::vector<double> const & mu,
                                       std::vector<double> const & y,
                                       std::vector<double> *fr) {


        //  for (std::ptrdiff_t i = 0; i < n0; ++i)
        //      for (std::size_t j = 0; j < nb; ++j)
        //          Bb.at(i*nb + j) = mu.at(n0-i+j);

        // cross_y[i] = produto interno da linha i de Bb com y descolado de jb+1
        // a linha i de Bb é dada por mu[n0-i+j] com j variando de 0 a nb-1
        // cross_y[i] = cblas_ddot(nb, mu.data() + n0 - i, 1,
        //                              y.data() + jb + 1, 1);

        for (std::size_t i = 0; i < n0; ++i) {
            (*fr)[i] = cblas_ddot(
                            nb,
                            mu.data() + n0 - i, // - 1,  // + n0 - i,
                            1,      // incremento em mu = 1
                            y.data() + jb + 1,
                            1);
        }
    }
#endif
