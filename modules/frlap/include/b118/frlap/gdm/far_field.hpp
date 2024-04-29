// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <utility>

namespace b118 { namespace frlap { namespace gdm {
namespace far_field {

    struct zero final {};

    template<class F>
    struct general final {
        F const y;
        general(F y) : y(y) {}
    };

    template<typename Real>
    struct algebraic final {
        Real const decay_minus;
        Real const decay_plus;
        Real const y_minus;
        Real const y_plus;
        Real const y_a;
        Real const y_b;
 
        algebraic(std::pair<Real, Real> decays,
                  std::pair<Real, Real> yinfty,
                  std::pair<Real, Real> ybndry)
            : decay_minus(decays.first), decay_plus(decays.second),
              y_minus    (yinfty.first), y_plus    (yinfty.second),
              y_a        (ybndry.first), y_b       (ybndry.second)
        {}

        algebraic(std::pair<Real, Real> decays,
                  Real yinfty,
                  std::pair<Real, Real> ybndry)
            : algebraic(decays, std::make_pair(yinfty, yinfty), ybndry)
        {}

        algebraic(std::pair<Real, Real> decays,
                  std::pair<Real, Real> ybndry)
            : algebraic(decays, static_cast<Real>(0), ybndry)
        {}

        algebraic(Real decay,
                  std::pair<Real, Real> yinfty,
                  std::pair<Real, Real> ybndry)
            : algebraic(std::make_pair(decay, decay), yinfty, ybndry)
        {
            std::cout << "Calling this constuctor" << std:: endl;
        }

        algebraic(Real decay,
                  std::pair<Real, Real> ybndry)
            : algebraic(decay, std::make_pair(static_cast<Real>(0),
                                              static_cast<Real>(0)), ybndry)
        {}
        
        algebraic(Real decay,
                  Real y_a, Real y_b)
            : algebraic(decay, std::make_pair(y_a, y_b))
        {}
    };

}    // end namespace far_field
}}}  // end namespace b118::frlap::gdm
