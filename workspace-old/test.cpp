// Copyright 2014 Guilherme F. Fornel <gffrnl@gmail.com>
// Compilar com:
// g++ -std=c++17 test.cpp -o test -Wall -Wpedantic -lfftw3f -lfftw3 -lfftw3l \
-lfftw3q -lopenblas -DNDEBUG

#include <iostream>
#include <cmath>
#include <b118/linspace.hpp>
#include <b118/closest_sorted.hpp>
#include <b118/frlap/gdm/coefficients/centered_3_point_periodized.hpp>
#include "./generalized_differences.hpp"

#define per3point centered_3_point_periodized

double const ealpha = 0.4;

double y(double x) {
    double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
    return sgn * std::pow(1.0 + x*x, -0.5 + ealpha/2.0);
}

double exact(double x) {
    double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
    return sgn * std::exp2(ealpha) *
        std::tgamma(0.5+ealpha/2.0) / std::tgamma(0.5-ealpha/2.0) *
        std::pow(1.0 + x*x, -0.5-ealpha/2.0);
}

int main() {
    using std::cout;
    using std::endl;
    using b118::exclude;
    namespace frlap = b118::frlap;
    namespace coeff = b118::frlap::gdm::coefficients;
    namespace far_field = b118::frlap::gdm::far_field;

    typedef frlap::generalized_differences<double, coeff::per3point> Method;

    double a  = -6.0;
    double b  = +6.0;
    double a0 = -2.0;
    double b0 = +2.0;
    std::size_t n = 121;
    double deltax;

    std::vector<double> Y(n);


    auto X =
        b118::linspace<std::vector, double>(a, b, n, exclude::none, &deltax);
    std::size_t ja = b118::closest_sorted(X.cbegin(), X.cend(),
                                          a0);
    std::size_t jb = b118::closest_sorted(X.cbegin(), X.cend(),
                                          b0);

    for (std::size_t k{0}; k < n; ++k) {
        Y.at(k) = y(X.at(k));
    }

    std::size_t n0 = jb - ja + 1;
    std::vector<double> FLY_0(n0);

    Method method(ealpha, a, b, n);
    method.compute_truncated(Y, ja, jb, FLY_0);

    std::vector<double> FF(FLY_0.size());
    frlap::gdm::far_field_estimator<double, far_field::general> ffest(
        y, ealpha, deltax, a, b);
    for (std::size_t j = ja; j <= jb; ++j) {
        FF[j-ja] = ffest(X[j]);
    }


    // Write the results to a file:
    if (true) {
        std::FILE *fp;

        if ((fp = std::fopen("bench-701-c3p.dat", "w")) == NULL) {
            (void) std::fprintf(stderr,
                                "error: fopen() :"
                                "cannot open file bench701.dat to write\n");
            std::abort();
        }

        (void) std::fprintf(fp,
                            "%-22s\t%-22s\t%-22s\t%-22s\t%-22s\t%-22s\t%-22s\n",
                            "X", "Y", "FLY0", "FF", "FLY", "Exact", "Error");
        for (std::size_t j = 0; j < n0; ++j)
            (void) std::fprintf(fp,
                                "%+22.15E\t%+22.15E\t%+22.15E\t%+22.15E\t"
                                "%+22.15E\t%+22.15E\t%+22.15E\n",
                                X[j+ja], Y[j+ja], FLY_0[j], FF[j],
                                FLY_0[j] + FF[j],
                                exact(X[j+ja]),
                                std::fabs(exact(X[j+ja])-(FLY_0[j] + FF[j])));
        std::fclose(fp);
    }

    return 0;
}
