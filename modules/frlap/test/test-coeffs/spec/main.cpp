// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <iostream>
#include "../../b118/frlap/generalized_differences_methods/coefficients/spectral.hpp"

namespace coefficients = b118::frlap::gdm::coefficients;

int main() {
    using std::cout;
    using std::endl;

    double ealpha = 1.4;
    double deltax = 0.1;
    std::size_t n = 100;

    coefficients::spectral<double> cgtor(ealpha, deltax);

    for (std::size_t k{0}; k < n; ++k)
        cout << "k : " << k << ", value: " << cgtor(k) << '\n';
    cout << endl;

    return 0;
}