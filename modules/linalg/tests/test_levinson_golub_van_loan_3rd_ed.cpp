// Copyright 2024 Guilherme F. Fornel

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "../include/b118/linalg/toeplitz/solvers/levinson.hpp"

int main() {
    using std::cout;
    using std::endl;
    using vec = std::vector<double>;
    auto print = [](double const & x) { cout << x << ' '; };

    std::cout.precision(8);

    // vec row = {1.0,  0.5, 0.2};
    // vec rhs = {4.0, -1.0, 3.0};
    // vec row = {2.0,  1.0, 0.4};
    // vec rhs = {8.0, -2.0, 6.0};
    vec row = {5.0,  2.5, 1.0};
    vec rhs = {20.0, -5.0, 15.0};
    vec sol_exact = {355.0/56, -376.0/56, 285.0/56};

    {
        cout << "sol_exact = ";
        std::for_each(sol_exact.begin(), sol_exact.end(), print);
        cout << endl;
    }

    vec sol(rhs.size());

    // -------------------------------------------------------------------------

    (void) b118::linalg::toeplitz::solvers::levinson<double>(row.begin(), row.end(), rhs.begin(), sol.begin());  // NOLINT
    cout << "sol       = ";
    std::for_each(sol.begin(), sol.end(), print);
    cout << endl;

    // std::reverse(row.begin(), row.end());
    // (void) levinson<double>(row.rbegin(), row.rend(), rhs.begin(), sol.begin());  // NOLINT 
    // cout << "sol       = ";
    // std::for_each(sol.begin(), sol.end(), print);
    // cout << endl;

    // -------------------------------------------------------------------------

    // cout << "sol       = ";
    // std::for_each(
    //     sol.begin(),
    //     levinson<double>(row.begin(), row.end(), rhs.begin(), sol.begin()),
    //     print);
    // cout << endl;

    // std::reverse(row.begin(), row.end());
    // cout << "sol       = ";
    // std::for_each(
    //     sol.begin(),
    //     levinson<double>(row.rbegin(), row.rend(), rhs.begin(), sol.begin()),
    //     print);
    // cout << endl;

    // -------------------------------------------------------------------------

    std::transform(sol_exact.cbegin(), sol_exact.cend(),
                   sol.cbegin(),
                   sol.begin(),
                   std::minus<double>{});
    std::copy(sol.begin(),
              sol.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

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
