// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>
#include <iostream>
#include "../../../modules/frlap/include/b118/frlap/gdm/coefficients/spectral_qawo.hpp"  // NOLINT
#include <format>

int main() {
    using std::cout;
    using std::endl;
    namespace coefficients = b118::frlap::gdm::coefficients;

    double ealpha  = 1.4;
    double deltax = 1;

    coefficients::spectral_qawo<double> cgtor(1000);
    cgtor.generate(ealpha, deltax);

    for (std::size_t k = 0; k < 1000; k+=11) {
        //cout << "mu[" << k << "] = " << cgtor.coeffs[k] << '\n';
        cout << std::format("mu[{}] = {:>20.17e}\n", k, cgtor.coeffs[k]);
    }
}