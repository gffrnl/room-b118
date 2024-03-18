// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <iostream>
#include <cmath>
#include "./grid.hpp"

int main() {
    using std::cout;
    using std::endl;

    grid<double> G1({-1.2, 2.1}, 12);
    grid<double> G2({-1.2, 2.1}, 12);

    cout << "G1 == G2 ? " << std::boolalpha << (G1 == G2) << '\n' << endl;

    cout << "G1.numnodes() = " << G1.numnodes() << endl;
    cout << "G1.spacing()  = " << G1.spacing()  << endl;

    cout << "G1.closest(-0.61) = " << G1.closest(-0.61) << endl;
    cout << "G1.closest(-0.46) = " << G1.closest(-0.46) << endl;
    cout << "G1.closest(-0.45) = " << G1.closest(-0.45) << endl;
    cout << "G1.closest(-0.44) = " << G1.closest(-0.44) << endl;

    for (grid<double>::size_type k{0}; k < 12; ++k)
        cout << G1[k] << ' ';
    cout << endl;

    double const ealpha = 0.4;
    auto y       = [&ealpha](double const & x) -> double {
        double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
        return sgn * pow(1.0 + x*x, -0.5 + ealpha/2.0);
    };

    grid<double> G3({-1.0, 1.0}, 11);
    grid_function<double> Y(G3, y);

    cout << "G3.numnodes() = " << G3.numnodes() << '\n';
    cout << "G3.spacing()  = " << G3.spacing() << '\n';

    // for (grid<double>::size_type k{0}; k < G3.numnodes(); ++k)
    //     cout << Y[k] << ' ';
    // cout << endl;

    for (grid<double>::size_type k{0}; k < G3.numnodes(); ++k)
        cout << G3[k] << ' ';
    cout << endl;

    grid<double>::size_type const ka = G3.closest(-0.61);
    grid<double>::size_type const kb = G3.closest(0.22);

    cout << "ka = " << ka << '\n';
    cout << "kb = " << kb << '\n';

    auto G4 = G3.subgrid(ka, kb);

    for (grid<double>::size_type k{0}; k < G4.numnodes(); ++k)
        cout << G4[k] << ' ';
    cout << endl;

    cout << "G4.numnodes() = " << G4.numnodes() << '\n';
    cout << "G4.spacing()  = " << G4.spacing() << '\n';

    cout << "has G3 G4 as subgrid? " << std::boolalpha << G3.has_subgrid(G4)
         << "\n";
    cout << "is G4 subgrid of G3? " << std::boolalpha << G4.is_subgrid(G3)
         << "\n";

    return 0;
}
