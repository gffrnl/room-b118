// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <b118/grid_function.hpp>
#include <b118/signal/convolution.hpp>
#include <b118/linalg/toeplitz/fast_symm_toeplitz_prod.hpp>
#include "./gdm/coefficients/generator.hpp"
#include "./gdm/coefficients/cper3point.hpp"
#include "./gdm/far_field.hpp"
#include "./gdm/far_field_estimator.hpp"

#include <iostream>
#include <cassert>

// // Convolution via dot product
// // The output fr satisfies:
// // fr[i] += sum(kernel[input_size - 1 - i + j] * y[j], j = 0 ... input_size - 1), i = 0 .. output_size - 1
// // kernel[k] must be defined in [0, output_size + input_size - 1]
// void conv1d_ingenua(std::size_t output_size, std::size_t input_size,
//                     Real const * const kernel,
//                     Real const * const y,
//                     Real       * const fr) {
//     // Obs: Blas shifts the pointer when INCX < 0:
//     //  ...
//     //     IX = 1
//     //     IY = 1
//     //     IF (INCX.LT.0) IX = (-N+1)*INCX + 1
//     //     IF (INCY.LT.0) IY = (-N+1)*INCY + 1
//     //     DO I = 1,N
//     //         DTEMP = DTEMP + DX(IX)*DY(IY)
//     //         IX = IX + INCX
//     //         IY = IY + INCY
//     //     END DO
//     // ...
//     for (std::size_t i = 0; i < output_size; ++i) {
//         fr[i] += cblas_ddot(input_size,
//                             const_cast<Real*>(kernel) + i,
//                             -1,      // increment on kernel = -1
//                             const_cast<Real*>(y),
//                             1);
//     }
// }


// // Cross-correlation via dot product
// // The output fr satisfies:
// // fr[i] += sum(kernel[output_size - 1 - i + j] * y[j], j = 0 ... input_size - 1), i = 0 .. output_size - 1
// // kernel[k] must be defined in [0, output_size + input_size - 1]
// void cross1d_ingenua(std::size_t output_size, std::size_t input_size,
//                      Real const * const kernel,
//                      Real const * const y,
//                      Real       * const fr) {
//     for (std::size_t i = 0; i < output_size; ++i) {
//         fr[i] += cblas_ddot(input_size,
//                             const_cast<Real*>(kernel) + (output_size - 1) - i,
//                             1,      // increment on kernel = 1
//                             const_cast<Real*>(y),
//                             1);
//     }
// }

namespace b118  { namespace frlap {

template<
    typename Real,
    template<typename> class CoeffGenerator =
        b118::frlap::gdm::coefficients::cper3point
    >
class generalized_differences final {
 public:
    generalized_differences(Real ealpha, b118::grid<Real> nodes)
        : m_ealpha(ealpha), G_(nodes)
    {
        m_deltax = G_[1]-G_[0];
    }
 
    generalized_differences(Real ealpha, Real a, Real b, std::size_t n)
        : m_ealpha(ealpha),
          m_deltax((b - a) / static_cast<Real>(n - 1))
    {}

    generalized_differences(Real ealpha, std::pair<Real, Real> ab,
        std::size_t n)
        : generalized_differences(ealpha, ab.first, ab.second, n)
    {}

    Real ealpha() const { return m_ealpha; }
    Real deltax() const { return m_deltax; }
    
    inline CoeffGenerator<Real> coefficients_generator() const {
        return CoeffGenerator<Real>(m_ealpha, m_deltax);
    }

    void compute_truncated(b118::grid_function<Real> const & Y,
                           b118::grid_function<Real>       * FLY0) {
        if (!((*FLY0).get_grid().is_subgrid(Y.get_grid())))
            throw std::invalid_argument("generalized_differences::"
                                        "compute_truncated() : "
                                        "the grid of *FLY0 is not a subgrid of"
                                        "Y");
        Real const xa = FLY0->get_grid()[0];
        Real const xb = FLY0->get_grid()[FLY0->get_grid().numnodes() - 1];
        std::size_t const ja = Y.get_grid().closest(xa);
        std::size_t const jb = Y.get_grid().closest(xb);

        // std::cout << "xa = " << xa << std::endl;
        // std::cout << "xb = " << xb << std::endl;
        // std::cout << "ja = " << ja << std::endl;
        // std::cout << "jb = " << jb << std::endl;

        std::size_t const n = Y.get_grid().numnodes();
        std::size_t const na = ja;
        std::size_t const nb = n-1-jb;
        std::size_t const n0 = jb-ja+1;
        std::size_t const nc = (jb+1 > n-ja)? jb+1 : n-ja;

        CoeffGenerator<Real> coeff_gtor = coefficients_generator();
        std::vector<Real> mu(nc);
        for (std::size_t k = 0; k < nc; ++k)
            mu.at(k) = coeff_gtor(k);
        
        // {
        //     std::size_t k = 0;
        //     std::for_each(mu.begin(), mu.end(),
        //                   [&coeff_gtor, &k](Real & x) -> Real {
        //                       x[k] = coeff_gtor(k++);
        //                   });
        // }

        // TODO(Fabio/Guilherme) criar condição para escolher entre convolução
        //                       e produto.
        // if (false) {
        //     conv1d_ingenua(n0, na,  m_cgtor.coeffs.data() + ja - na + 1, y.data(),
        //                             frLap_y.data());

        //     cross1d_ingenua(n0, nb, m_cgtor.coeffs.data() + 1, y.data() + jb + 1,
        //                             frLap_y.data());
        // } else {
        std::size_t const conv_size = n0 + std::max(na, nb) - 1;
        b118::signal::convolution<Real> conv(conv_size);
        conv.create_plans(conv_size);
        conv.conv(n0, na, mu.data() + ja - na + 1, Y.data(),
                  FLY0->data(),
                  false);
        conv.conv(n0, nb, mu.data() + 1, Y.data() + jb + 1,
                  FLY0->data(),
                  true);
        // }


        // 5.4. Fast symmetric toeplitz-vector A*Yint:
        b118::linalg::toeplitz::fast_symm_toeplitz_prod<Real>{}(
            mu.data(), mu.data() + n0,
                                         Y.data()+ja, FLY0->data(),
                                         true);

    }


    // inline void compute_truncated(b118::grid_function<Real> const & Y,
    //                               b118::grid_function<Real>       * FLY0) {
    //     compute_truncated(coefficients_generator(), Y, FLY0);
    // }
        
    // void compute_truncated(b118::grid_function<Real> const & Y,
    //                        b118::grid_function<Real>       * FLY0) {
    //     if (!((*FLY0).get_grid().is_subgrid(Y.get_grid())))
    //         throw std::invalid_argument("generalized_differences::"
    //                                     "compute_truncated() : "
    //                                     "the grid of *FLY0 is not a subgrid of"
    //                                     "Y");
    //     Real const xa = FLY0->get_grid()[0];
    //     Real const xb = FLY0->get_grid()[FLY0->get_grid().numnodes() - 1];
    //     std::size_t const ja = Y.get_grid().closest(xa);
    //     std::size_t const jb = Y.get_grid().closest(xb);

    //     // std::cout << "xa = " << xa << std::endl;
    //     // std::cout << "xb = " << xb << std::endl;
    //     // std::cout << "ja = " << ja << std::endl;
    //     // std::cout << "jb = " << jb << std::endl;

    //     std::size_t const n = Y.get_grid().numnodes();
    //     std::size_t const na = ja;
    //     std::size_t const nb = n-1-jb;
    //     std::size_t const n0 = jb-ja+1;
    //     std::size_t const nc = (jb+1 > n-ja)? jb+1 : n-ja;

    //     auto coeff_gtor = coefficients_generator();
    //     std::vector<Real> mu(nc);
    //     for (std::size_t k = 0; k < nc; ++k)
    //         mu.at(k) = coeff_gtor(k);
        
    //     // {
    //     //     std::size_t k = 0;
    //     //     std::for_each(mu.begin(), mu.end(),
    //     //                   [&coeff_gtor, &k](Real & x) -> Real {
    //     //                       x[k] = coeff_gtor(k++);
    //     //                   });
    //     // }

    //     // TODO(Fabio/Guilherme) criar condição para escolher entre convolução
    //     //                       e produto.
    //     // if (false) {
    //     //     conv1d_ingenua(n0, na,  m_cgtor.coeffs.data() + ja - na + 1, y.data(),
    //     //                             frLap_y.data());

    //     //     cross1d_ingenua(n0, nb, m_cgtor.coeffs.data() + 1, y.data() + jb + 1,
    //     //                             frLap_y.data());
    //     // } else {
    //     std::size_t const conv_size = n0 + std::max(na, nb) - 1;
    //     convolution<Real> conv(conv_size);
    //     conv.create_plans(conv_size);
    //     conv.conv(n0, na, mu.data() + ja - na + 1, Y.data(),
    //               FLY0->data(),
    //               false);
    //     conv.conv(n0, nb, mu.data() + 1, Y.data() + jb + 1,
    //               FLY0->data(),
    //               true);
    //     // }


    //     // 5.4. Fast symmetric toeplitz-vector A*Yint:
    //     fast_symm_toeplitz_product<Real>(mu.data(), mu.data() + n0,
    //                                      Y.data()+ja, FLY0->data(),
    //                                      true);

    // }


    template<
        template<class...> class FarFieldEstimatorKind,
        class ...Args
        >
    void compute_far_field(FarFieldEstimatorKind<Args...>   ffkind,
                           b118::grid_function<Real>      * FF    ,
                           bool inplace = false) {
        b118::frlap::gdm::far_field_estimator<
            Real, FarFieldEstimatorKind
            > ffest(m_ealpha, ffkind, G_, FF->get_grid());
        auto const the_grid = FF->get_grid();
        // TODO(gffrnl): check the_grid == G0_;
        assert(the_grid.is_subgrid(G_));
        auto const numnodes = the_grid.numnodes();
        if (!inplace) {
            for (std::size_t k = 0; k < numnodes; ++k)
                (*FF)[k] = ffest(the_grid[k]);
        } else {
            for (std::size_t k = 0; k < numnodes; ++k)
                (*FF)[k] += ffest(the_grid[k]);
        }
    }

    template<
        template<class...> class FarFieldEstimatorKind,
        class ...Args
        >
    void compute(FarFieldEstimatorKind<Args...>         ffkind,
                 b118::grid_function<Real>      const & Y,
                 b118::grid_function<Real>            * FL    ) {
        compute_truncated(Y, FL);
        compute_far_field(ffkind, FL, true);
    }
    
 private:
    Real m_ealpha;
    Real m_deltax;

    b118::grid<Real> G_;
    b118::grid<Real> G0_;
};


}}  // end namespace b118::frlap

