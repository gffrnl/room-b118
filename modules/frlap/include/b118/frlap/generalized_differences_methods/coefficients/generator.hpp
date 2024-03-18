// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>
// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once

#include <cstddef>
#include <vector>

namespace b118         {
namespace frlap        {
namespace gdm          {
namespace coefficients {

template<typename Real, template<typename> class Method>
class generator {
 public:
    using method = Method<Real>;

    generator(Real ealpha, Real deltax)
        : ealpha(ealpha), deltax(deltax)
    {}

    Real operator() (std::size_t k) {
        return reinterpret_cast<Method<Real> *>(this)->operator()(k);
    }

 protected:
    Real ealpha;
    Real deltax;
};

}  // end namespace coefficients
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
