/*   libb118
 *
 *   modules/frlap/include/b118/frlap/gdm/far_field_estimator.hpp
 *   
 *   Truncated uniform gdm
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

#pragma once

#include <limits>
#include <functional>
#include <b118/frlap.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace b118 {
namespace frlap {
namespace gdm {

struct general {};
struct algebraic_decay {};

template<class EstimatorType>
struct far_field_estimator;

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118



namespace b118 {
namespace frlap {
namespace gdm {

template<>
struct far_field_estimator<general> {
    struct integrand_minus final {
        double ealpha;
        std::function<double(double)> y;
     public:
        double xj;

        template<class F>
        integrand_minus(double ealpha, F y)
            : ealpha(ealpha), y(y), xj(0)
        {}

        double operator()(double x) const {
            return y(x) / std::pow(xj - x, ealpha+1);
        }
    };

    struct integrand_plus final {
        double ealpha;
        std::function<double(double)> y;
     public:
        double xj;

        template<class F>
        integrand_plus(double ealpha, F y)
            : ealpha(ealpha), y(y), xj(0)
        {}

        double operator()(double x) const {
            return y(x) / std::pow(x - xj, ealpha+1);
        }
    };

    double ealpha;
    double limm, limp;
    double C1a;

    integrand_minus fm;
    integrand_plus  fp;

    boost::math::quadrature::tanh_sinh<double> integrator;

    template<class F>
    far_field_estimator(F y, double ealpha,
                        double limm = 0,
                        double limp = 0)
        : limm(limm),
          limp(limp),
          C1a(b118::frlap::normal_const<1>(ealpha)),
          fm(integrand_minus(ealpha, y)),
          fp(integrand_plus(ealpha, y))
    {}

    double operator()(double x) {
        fm.xj = x;
        fp.xj = x;
        return - C1a * (
            integrator.integrate(fm,
              - std::numeric_limits<double>::infinity(),
                limm)
            +
            integrator.integrate(fp,
                limp,
                std::numeric_limits<double>::infinity()));
    }
};

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118


namespace b118 {
namespace frlap {
namespace gdm {


template<>
struct far_field_estimator<algebraic_decay> {

    double ealpha;
    double edecay;
    double xa;
    double xb;
    double xja;
    double xjb;
    double deltax;
    double yja;
    double yjb;
    double C1a;

    template<class F>
    far_field_estimator(F y,
                        double ealpha,
                        double edecay,
                        double xa,
                        double xb,
                        double xja,
                        double xjb,
                        double deltax)
        : ealpha(ealpha),
          edecay(edecay),
          xa(xa),
          xb(xb),
          xja(xja),
          xjb(xjb),
          deltax(deltax),
          yja(y(xja)),
          yjb(y(xjb)),
          C1a(b118::frlap::normal_const<1>(ealpha))
    {}

    double operator()(double x) {
        using boost::math::hypergeometric_pFq;
        return -C1a/(ealpha+edecay) * (
            yja * std::pow(std::fabs(xja), edecay)
                 / std::pow(std::fabs(xa - deltax/2), ealpha + edecay)
                 * hypergeometric_pFq({ealpha + 1, ealpha + edecay},
                                      {ealpha + edecay + 1},
                                      - x / std::fabs(xa - deltax/2))
            +
            yjb * std::pow(std::fabs(xjb), edecay)
                 / std::pow(std::fabs(xb + deltax/2), ealpha + edecay)
                 * hypergeometric_pFq({ealpha + 1, ealpha + edecay},
                                      {ealpha + edecay + 1},
                                      + x / std::fabs(xb + deltax/2)) );
    }
};

}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118

