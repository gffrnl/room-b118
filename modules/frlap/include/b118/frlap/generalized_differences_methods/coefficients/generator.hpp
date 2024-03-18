// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once

#include <cstddef>
#include <vector>

namespace b118         {
namespace frlap        {
namespace gdm          {
namespace coefficients {

template<typename Real, template<typename> class Method>
struct generator {
    using method = Method<Real>;

    std::vector<Real> coeffs;

    explicit generator(std::size_t n)
        : coeffs(std::vector<Real>(n))
    {}

    void generate(Real ealpha, Real deltax) {
        reinterpret_cast<Method<Real> *>(this)->generate(ealpha, deltax);
    }
};

}  // end namespace coefficients
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
