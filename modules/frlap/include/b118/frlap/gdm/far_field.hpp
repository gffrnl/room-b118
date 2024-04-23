// Copyright 2024 Guilherme F. Fornel

#pragma once

namespace b118 { namespace frlap { namespace gdm {

/*
template<class Derived>
struct far_field_kind {
  typedef Derived kind;
};
*/

namespace far_field {

struct zero final {}; //: public b118::frlap::gdm::far_field_kind<far_field::zero> {};

template<class F>
struct general final { //: public far_field_kind<far_field::general>  {
    F const y;
    general(F y) : y(y) {}
};

template<typename Real>
struct algebraic final {//: public far_field_kind<far_field::algebraic> {
    Real const edecay;
    Real const yja, yjb;
    algebraic(Real edecay, Real yja, Real yjb)
        : edecay(edecay), yja(yja), yjb(yjb)
    {}
};

}    // end namespace far_field
}}}  // end namespace b118::frlap::gdm
