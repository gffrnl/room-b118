// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>
// Compile with:
// g++ -g -std=c++17 -Wall -Wpedantic main.cpp -o b118-linalg-fsymmtp-example
//    -lfftw3f -lfftw3 -lfftw3l -lfftw3q

#include <iostream>
#include <vector>
#include <algorithm>
#include "../include/b118/linalg/toeplitz/fast_symm_toeplitz_product.hpp"

template<typename Real>
auto print = [](Real const & x) {
    std::cout << x << ' ';
};

template<>
auto print<__float128> = [](__float128 const & x) {
    std::cout << static_cast<long double>(x) << ' ';
};

int main() {
    using std::cout;
    using std::clog;
    using std::endl;
    using real_t = __float128;

    clog << "b118::linalg :"
         << "fast symmetric Toeplitz product example"
         << endl;

    std::vector<real_t>       row1 = {5.0,  2.5, 1.0};
    std::vector<real_t>       x    = {355.0/56, -376.0/56, 285.0/56};
    std::vector<real_t> const b    = {20.0, -5.0, 15.0};
    std::vector<real_t>       ba(x.size());


    cout << "b  = ";
    std::for_each(b.begin(), b.end(), print<real_t>);
    cout << endl;

    cout << "ba = ";
    std::for_each(ba.begin(), ba.end(), print<real_t>);
    cout << endl;

    fast_symm_toeplitz_product<real_t>(row1.cbegin(), row1.cend(),
                                       x.cbegin(),
                                       ba.begin());

    cout << "ba = ";
    std::for_each(ba.begin(), ba.end(), print<real_t>);
    cout << endl;

    return 0;
}
