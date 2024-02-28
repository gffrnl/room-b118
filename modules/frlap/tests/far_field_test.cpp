// Copyright 2014 Guilherm F. Fornel <gffrnl@gmail.com>

#include <iostream>
#include <b118/frlap/gdm/far_field_estimator.hpp>

int main() {
    using std::cout;
    using std::endl;
    using namespace b118::frlap::gdm;  // NOLINT

    cout << "Example" << endl;

    double ealpha, edecay;
    ealpha  = 1.4;
    edecay = 0.8;
    double xj = 0.0;
    auto y = [ealpha, edecay](double x) -> double {
        return std::pow(std::abs(x), -(ealpha+edecay));
    };

    far_field_estimator<general> ffest(y, 0.4, -6.0, 6.0);

    cout << "ffest(0.0)  = " << ffest(0.0)  << endl;

    cout << "ffest(-2.0) = " << ffest(-2.0) << endl;
    cout << "ffest(+2.0) = " << ffest(+2.0) << endl;

    cout << "ffest(-5.0) = " << ffest(-5.9) << endl;
    cout << "ffest(+5.0) = " << ffest(+5.9) << endl;

    cout << "-----------------------------------------" << endl;

    far_field_estimator<algebraic_decay> ffest2(y, 0.4, ealpha+edecay,
                                                -6.0, 6.0,
                                                -2.0, 2.0,
                                                0.1);

    cout << "ffest2(0.0)  = " << ffest2(0.0)  << endl;

    cout << "ffest2(-2.0) = " << ffest2(-2.0) << endl;
    cout << "ffest2(+2.0) = " << ffest2(+2.0) << endl;

    cout << "ffest2(-5.0) = " << ffest2(-5.9) << endl;
    cout << "ffest2(+5.0) = " << ffest2(+5.9) << endl;
    return 0;
}


// Expected result:
//
// ffest(0.0)  = -0.00121056
// ffest(-2.0) = -0.00135196
// ffest(+2.0) = -0.00135196
// ffest(-5.0) = -0.0138671
// ffest(+5.0) = -0.0138671
// -----------------------------------------
// ffest2(0.0)  = -0.00118472
// ffest2(-2.0) = -0.00132058
// ffest2(+2.0) = -0.00132058
// ffest2(-5.0) = -0.0108574
// ffest2(+5.0) = -0.0108574
