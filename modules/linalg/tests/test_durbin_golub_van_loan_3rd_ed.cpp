// Copyright 2024 Guilherme F. Fornel

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "../include/b118/linalg/toeplitz/solvers/durbin.hpp"

int main() {
    using vec = std::vector<double>;
    auto print = [](double const & x) { std::cout << x << ' '; };

    std::cout.precision(8);

    vec rhs = {0.5, 0.2, 0.1};
    vec sol_exact = {-75.0/140, 12.0/140, -5.0/140};
    std::cout << "sol_exact = ";
    std::for_each(sol_exact.begin(), sol_exact.end(), print);
    std::cout << std::endl;


    vec sol(rhs.size());

    std::cout << "sol       = ";
    b118::linalg::toeplitz::solvers::durbin<double>(rhs.begin(), rhs.end(),
                                            sol.begin());
    std::for_each(sol.begin(), sol.end(), print);
    // std::for_each(sol.begin(),
    //               b118::linalg::toeplitz::solvers::durbin<double>(rhs.begin(),
    //                                                       rhs.end(),
    //                                                       sol.begin()),
    //               print);
    std::cout << std::endl;

    std::transform(sol_exact.cbegin(), sol_exact.cend(),
                   sol.cbegin(),
                   sol.begin(),
                   std::minus<double>{});
    // std::copy(sol.begin(),
    //           sol.end(),
    //           std::ostream_iterator<double>(std::cout, " "));
    // std::cout << std::endl;

    double maxabserr = 1;
    auto [min, max] =
        std::minmax_element(sol.begin(), sol.end());
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
