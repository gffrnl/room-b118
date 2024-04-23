// Copyright 2024 Guilherme F. Fornel

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "../include/b118/linalg/toeplitz/solvers/trench.hpp"

int main() {
    namespace toeplitz = b118::linalg::toeplitz;
    using vec = std::vector<double>;
    auto print = [](double const & x) { std::cout << x << ' '; };

    std::cout.precision(8);

    vec row = {1, 0.5, 0.2};
    vec inv_exact = {75.0/56, -5.0/7, 5.0/56, 12.0/7};
    std::cout << "inv_exact = ";
    std::for_each(inv_exact.begin(), inv_exact.end(), print);
    std::cout << std::endl;


    vec inv(inv_exact.size());

    std::cout << "inv       = ";
    toeplitz::solvers::trench<double>(row.begin(), row.end(), inv.begin());
    std::for_each(inv.begin(), inv.end(), print);
    // std::for_each(inv.begin(),
    //               toeplitz::solvers::trench<double>(row.begin(), row.end(),
    //                                                 inv.begin()),
    //               rint);
    std::cout << std::endl;

    std::transform(inv_exact.cbegin(), inv_exact.cend(),
                   inv.cbegin(),
                   inv.begin(),
                   std::minus<double>{});
    // std::copy(inv.begin(),
    //           inv.end(),
    //           std::ostream_iterator<double>(std::cout, " "));
    // std::cout << std::endl;

    double maxabserr = 1;
    auto [min, max] =
        std::minmax_element(inv.begin(), inv.end());
    std::cout << "*min      = " << *min << std::endl;
    std::cout << "*max      = " << *max << std::endl;
    if (std::fabs(*min) > std::fabs(*max))
        maxabserr = std::fabs(*min);
    else
        maxabserr = std::fabs(*max);
    std::cout << "maxabserr = " << maxabserr << std::endl;

    if (maxabserr > 1.0e-6)
        return 1;  // Test not passed

    return 0;  // Test passed!
}
