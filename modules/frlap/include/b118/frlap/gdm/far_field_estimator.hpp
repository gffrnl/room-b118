// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <functional>
#include <b118/grid.hpp>
#include "../../frlap.hpp"
#include "./far_field.hpp"
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace b118 { namespace frlap {
namespace gdm {

    template<typename Real, template<class...> class EstimatorKind>
    class far_field_estimator;


    // template<typename Real>
    // class far_field_estimator<Real, far_field::zero> final {
    //  public:
    //     far_field_estimator() {}
    //     Real operator()([[maybe_unused]] Real const & x) {
    //         return static_cast<Real>(0);
    //     }
    // };

    
    template<typename Real>
    class far_field_estimator<Real, far_field::general> final {
     public:
        template<class F>
        far_field_estimator(Real ealpha,
                            far_field::general<F> ffkind,
                            b118::grid<Real> G)
            : m_deltax(G[1]- G[0]),
              m_a(G[0]),
              m_b(G[G.numnodes() - 1]),
              m_C1a(b118::frlap::normal_const<1>(ealpha)),
              m_fminus(integrand_minus(ealpha, m_deltax, ffkind.y)),
              m_fplus (integrand_plus (ealpha, m_deltax, ffkind.y))  // NOLINT
        {}

        Real operator()(Real x) {
            m_fminus.xj = x;
            m_fplus.xj  = x;
            return - m_C1a * (integrator.integrate(m_fminus,
                             - std::numeric_limits<Real>::infinity(),
                              m_a - m_deltax / static_cast<Real>(2))
            +             integrator.integrate(m_fplus,
                              m_b + m_deltax / static_cast<Real>(2),
                              std::numeric_limits<Real>::infinity()));
        }

     private:
        struct integrand_minus final {
            Real ealpha;
            Real deltax;
            std::function<Real(Real)> y;
         public:
            Real xj;

            template<class F>
            integrand_minus(Real ealpha, Real deltax, F y) : ealpha(ealpha),
                                                             deltax(deltax),
                                                             y(y),
                                                             xj(0)
            {}

            Real operator()(Real x) const {
                return y(x) / pow(xj - x, ealpha + static_cast<Real>(1));
            }
        };

        struct integrand_plus final {
            Real ealpha;
            Real deltax;
            std::function<Real(Real)> y;
         public:
            Real xj;

            template<class F>
            integrand_plus(Real ealpha, Real deltax, F y)  : ealpha(ealpha),
                                                             deltax(deltax),
                                                             y(y),
                                                             xj(0)
            {}

            Real operator()(Real x) const {
                return y(x) / std::pow(x - xj, ealpha + static_cast<Real>(1));
            }
        };

     private:
        Real m_ealpha;
        Real m_deltax;
        Real m_a, m_b;
        Real m_C1a;

        integrand_minus m_fminus;
        integrand_plus  m_fplus;

        boost::math::quadrature::tanh_sinh<Real> integrator;
    };


    template<typename Real>
    class far_field_estimator<Real, far_field::algebraic> final {
     public:
        far_field_estimator(Real alpha,
                            far_field::algebraic<Real> ffkind,
                            b118::grid<Real> GhJ)
            : alpha(alpha),
              invalpha(static_cast<Real>(1) / alpha),
              C_1alpha(b118::frlap::normal_const<1>(alpha)),
              c1(alpha + static_cast<Real>(1)),
              c2_minus(alpha + ffkind.decay_minus),
              c2_plus (alpha + ffkind.decay_plus ),
              c3_minus(alpha + ffkind.decay_minus + static_cast<Real>(1)),
              c3_plus (alpha + ffkind.decay_plus  + static_cast<Real>(1)),
              y_minus (ffkind.y_minus),
              y_plus  (ffkind.y_plus)
        {
            Real const h = GhJ.spacing();
            Real const a = GhJ.lendpoint();
            Real const b = GhJ.rendpoint();
            
            assert((a - h / static_cast<Real>(2)) < static_cast<Real>(0));
            assert((b + h / static_cast<Real>(2)) > static_cast<Real>(0));   
            
            xah2 = -a + h / static_cast<Real>(2);
            xbh2 =  b + h / static_cast<Real>(2);

            k_minus = ffkind.y_a * pow(-a, ffkind.decay_minus)
                    / ( c2_minus + pow(xah2, c2_minus) );
            k_plus  = ffkind.y_b * pow( b, ffkind.decay_plus )
                    / ( c2_plus  + pow(xbh2, c2_plus ) );

            /*
            std::cout << "alpha    = " << alpha << std::endl;
            std::cout << "invalpha = " << invalpha << std::endl;
            std::cout << "C_1alpha = " << C_1alpha << std::endl;
            std::cout << "c1       = " << c1 << std::endl;
            std::cout << "c2_minus = " << c2_minus << std::endl;
            std::cout << "c2_plus  = " << c2_plus << std::endl;
            std::cout << "c3_minus = " << c3_minus << std::endl;
            std::cout << "c3_plus  = " << c3_plus << std::endl;
            std::cout << "y_minus  = " << y_minus << std::endl;
            std::cout << "y_plus   = " << y_plus << std::endl;

            std::cout << "-----------------------" << std::endl;
            std::cout << "h        = " << h << std::endl;
            std::cout << "a        = " << a << std::endl;
            std::cout << "b        = " << b << std::endl;
            std::cout << "xah2     = " << xah2 << std::endl;
            std::cout << "xbh2     = " << xbh2 << std::endl;
            std::cout << "k_minus  = " << k_minus << std::endl;
            std::cout << "k_plus   = " << k_plus << std::endl;
            */
        }

        Real operator()(Real x) {
            return -C_1alpha * (
                (y_minus * Ih_minus(x) + k_minus * hyperg21_minus(x))
                +
                (y_plus  * Ih_plus (x) + k_plus  * hyperg21_plus (x))
            );
        }

     private:
        Real const alpha;
        Real const invalpha;  // 1 / alpha
        Real const C_1alpha;  // frac. Laplacian normalization constant
        Real const c1;        // alpha + 1
        Real const c2_minus;  // alpha + beta^{-}
        Real const c2_plus;   // alpha + beta^{+}
        Real const c3_minus;  // alpha + beta^{-} + 1
        Real const c3_plus;   // alpha + beta^{+}  + 1
        Real const y_minus;   // y^{-\infty}
        Real const y_plus;    // y^{+\infty}
        
        Real xah2;  // -a + h/2
        Real xbh2;  //  b + h/2

        Real k_minus;  // y(a) * (-a)^{beta^{-}} /
                       // ( (alpha +beta^{-}) + (-a + h/2)^{alpha + beta^{-}})
        Real k_plus;   // y(b) *    b^{beta^{+}} /
                       // ( (alpha +beta^{+}) + ( b + h/2)^{alpha + beta^{+}})

        inline Real Ih_minus(Real const & x) const {
            return invalpha * pow(xah2 + x, -alpha);
        }

        inline Real Ih_plus (Real const & x) const {
            return invalpha * pow(xbh2 - x, -alpha);
        }

        inline Real hyperg21_minus(Real const & x) const {
            using boost::math::hypergeometric_pFq;
            return hypergeometric_pFq({c1, c2_minus}, {c3_minus}, - x / xah2);
        }

        inline Real hyperg21_plus (Real const & x) const {
            using boost::math::hypergeometric_pFq;
            return hypergeometric_pFq({c1, c2_plus }, {c3_plus }, + x / xbh2); 
        }
    };

}  // end namespace gdm 
}} // end namespace b118::frlap
