// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <cmath>
#include <limits>
#include <functional>
#include <b118/grid.hpp>
#include "../../frlap.hpp"
#include "./far_field.hpp"
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace b118 { namespace frlap { namespace gdm {

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
                        b118::grid<Real> G,
                        b118::grid<Real> G0)
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
    far_field_estimator(Real ealpha,
                        far_field::algebraic<Real> ffkind,
                        b118::grid<Real> domain,
                        b118::grid<Real> inner_domain)
        : m_ealpha(ealpha),
          m_deltax(domain[1] - domain[0]),
          m_xa(domain[0]),
          m_xb(domain[domain.numnodes() - 1]),
          m_xja(inner_domain[0]),
          m_xjb(inner_domain[inner_domain.numnodes() - 1]),
          m_edecay(ffkind.edecay),
          m_yja(ffkind.yja),
          m_yjb(ffkind.yjb),
          m_C1a(b118::frlap::normal_const<1>(ealpha))
    {}

    Real operator()(Real x) {
        using boost::math::hypergeometric_pFq;
        return -m_C1a / (m_ealpha + m_edecay) * (
            m_yja * pow(fabs(m_xja), m_edecay)
                  / pow(fabs(m_xa - m_deltax / static_cast<Real>(2)),
                        m_ealpha + m_edecay)
                  * hypergeometric_pFq({m_ealpha + static_cast<Real>(1),
                                        m_ealpha + m_edecay},
                                       {m_ealpha + m_edecay
                                                 + static_cast<Real>(1)},
                                       - x / fabs(m_xa - m_deltax /
                                                        static_cast<Real>(2)))
            +
            m_yjb * pow(fabs(m_xjb), m_edecay)
                  / pow(fabs(m_xb + m_deltax / static_cast<Real>(2)),
                             m_ealpha + m_edecay)
                  * hypergeometric_pFq({m_ealpha + static_cast<Real>(1),
                                        m_ealpha + m_edecay},
                                       {m_ealpha + m_edecay
                                                 + static_cast<Real>(1)},
                                       + x / fabs(m_xb + m_deltax /
                                                        static_cast<Real>(2))));
    }

 private:
    Real m_ealpha;
    Real m_deltax;
    Real m_edecay;
    Real m_xa;
    Real m_xb;
    Real m_xja;
    Real m_xjb;
    Real m_yja;
    Real m_yjb;
    Real m_C1a;
};

}}} // end namespace b118::frlap::gdm
