// Copyright 2024 Guilherme F. Fornel

#pragma once

namespace b118      {
namespace frlap     {
namespace gdm       {
namespace far_field {

struct zero      {};
struct general   {};
struct algebraic {};

}  // end namespace far_field
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118

namespace b118 {
namespace frlap {
namespace gdm {

template<typename Real, class EstimatorType>
class far_field_estimator;

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118

#include <cmath>
#include <limits>
#include <functional>
#include <b118/frlap.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace b118 {
namespace frlap {
namespace gdm {

template<typename Real>
class far_field_estimator<Real, far_field::zero> final {
 public:
    // template<class F>
    // far_field_estimator(F y, Real ealpha, Real deltax, Real a, Real b) {}
    far_field_estimator() {}
    Real operator()(Real x) { return static_cast<Real>(0); }
};

template<typename Real>
class far_field_estimator<Real, far_field::general> final {
 public:
    template<class F>
    far_field_estimator(F y, Real ealpha, Real deltax, Real a, Real b)
        : m_deltax(deltax),
          m_a(a),
          m_b(b),
          m_C1a(b118::frlap::normal_const<1>(ealpha)),
          m_fminus(integrand_minus(ealpha, deltax, y)),
          m_fplus (integrand_plus (ealpha, deltax, y))  // NOLINT
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
    template<class F>
    far_field_estimator(F y,
                        Real ealpha, Real deltax, Real xa, Real xb,
                        Real xja, Real xjb, Real edecay)
        : m_ealpha(ealpha),
          m_deltax(deltax),
          m_xa(xa),
          m_xb(xb),
          m_xja(xja),
          m_xjb(xjb),
          m_edecay(edecay),
          m_yja(y(xja)),
          m_yjb(y(xjb)),
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

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
