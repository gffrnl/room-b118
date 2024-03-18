// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <iostream>
#include <cmath>
#include "../../grid.hpp"
#include "../../b118/frlap/generalized_differences.hpp"
#include "../../b118/linalg/vector.hpp"
#include "../../b118/linalg/matrix/symmetric_toeplitz.hpp"
#include "../../b118/linalg/toeplitz/fast_symm_toeplitz_product.hpp"
namespace linalg = b118::linalg;
namespace frlap  = b118::frlap;

int main() {
    using std::cout;
    using std::endl;
    using linalg::matrix_kind::symmetric_toeplitz;

    auto u0 = [](double const & x) -> double {
        return exp(-x*x);
    };

    double ealpha = 0.4;
    double kappa  = 1.0;

    double time = 0;
    double final_time = 10;
    double time_step = 0.001;

    double a = -1.0;
    double b = +1.0;
    double n =  5;

    frlap::generalized_differences<double> method(ealpha, {a, b}, n);

    // contructs the vector with initial conditions
    linalg::vector<double> U0(n);
    linalg::vector<double> U1(n);
    {
        grid<double> G({a, b}, n);
        // for (std::size_t k{0}; k < n; ++k)
        //     cout << G[k] << " .. ";
        cout << endl;
        for (std::size_t k{1}; k <= n; ++k)
            U0(k) = u0(G[k - 1]);
    }
    cout << "U0 = " << U0 << endl;

    // constructs the toeplitz matrix Hm1
    linalg::matrix<double, symmetric_toeplitz> Hm1(n);
    {
        double const mult = kappa * time_step;
        auto cgtor = method.get_cgtor();
        cgtor.coeffs.resize(n);
        cgtor.generate(ealpha, method.deltax());
        // cout << "coefficients:\n";
        // for (auto mu : cgtor.coeffs)
        //     cout << mu << ", ";
        // cout << endl;
        Hm1(1, 1) = 1.0 - mult * cgtor.coeffs.at(0);
        for (std::size_t k{2}; k <= n; ++k)
            Hm1(1, k) = - mult * cgtor.coeffs.at(k - 1);
    }

    while ((time += time_step) < final_time) {
        cout << "\ntime: " << time << endl;

        fast_symm_toeplitz_product<double>(Hm1.begin(), Hm1.end(),
                                           U0.begin(),
                                           U1.begin());
        cout << "U = " << U1 << endl;
        U0 = U1;
    }
    return 0;
}
