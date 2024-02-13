// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon (12-fev-2023)

#include "../equal_within_ulps.hpp"
#include <iostream>
#include <iomanip>

int main() {
    using namespace cppreference;
    double x = 0.3;
    double y = 0.1 + 0.2;
    std::cout << std::hexfloat;
    std::cout << "x = " << x << '\n';
    std::cout << "y = " << y << '\n';
    std::cout << (x == y ? "x == y" : "x != y") << '\n';
    for (std::size_t n = 0; n <= 10; ++n)
        if (equal_within_ulps(x, y, n)) {
            std::cout << "x equals y witshin " << n << " ulps" << '\n';
            break;
        }
}