// Copyright 2024 Guilherme F. Fornel

#pragma once

namespace b118 { namespace frlap { namespace gdm {
namespace far_field {

struct zero      {};

template<class F>
struct general final {
    F y;
    general(F y) : y(y) {}
};

template<typename Real>
struct algebraic {
    Real edecay;
    Real yja, yjb;
    algebraic(Real edecay, Real yja, Real yjb)
        : edecay(edecay), yja(yja), yjb(yjb)
    {}
};

}    // end namespace far_field
}}}  // end namespace b118::frlap::gdm
