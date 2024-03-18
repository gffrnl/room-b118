// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <cstddef>
#include <cmath>
#include <iostream>
#include <string>
#include "../install/include/b118/linspace.hpp"
#include "../install/include/b118/frlap/generalized_differences.hpp"
#include "../install/include/b118/frlap/generalized_differences_methods/coefficients/huang_oberman_2.hpp"
#include "../install/include/b118/frlap/generalized_differences_methods/far_field.hpp"

void write_results(std::string filename,
                   std::size_t n0,
                   std::size_t ja,
                   double const * const X,
                   double const * const Y,
                   double const * const FLY_0,
                   double const * const FF,
                   double const * const Exact);

namespace frlap        = b118::frlap;
namespace gdm          = b118::frlap::gdm;
namespace coefficients = b118::frlap::gdm::coefficients;
namespace far_field    = b118::frlap::gdm::far_field;

int main() {
    using std::cout;
    using std::endl;
    double const ealpha = 0.4;
    double const deltax = 0.1;
    double const a = -6.0;
    double const b = +6.0;
    double const a0 = -2.0;
    double const b0 = +2.0;
    double const n = 121;
    double const edecay = 0.8;

    auto y       = [&ealpha](double const & x) -> double {
        double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
        return sgn * pow(1.0 + x*x, -0.5 + ealpha/2.0);
    };
    auto frLap_y = [&ealpha](double const & x) -> double {
        double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
        return sgn * pow(1.0 + x*x, -0.5 - ealpha/2.0) * exp2(ealpha)
            * tgamma(0.5+ealpha/2.0) / tgamma(0.5-ealpha/2.0);
    };

    frlap::generalized_differences method(ealpha, {a, b}, n);
    std::vector<double> Y(n);
    b118::linspace(Y.begin(), Y.end(), a, b);
    std::vector<double> FLYe{Y};

    std::for_each(Y.begin(), Y.end(), [y](double & x) { x = y(x); });
    // std::for_each(FLYe.begin(), FLYe.end(), frLap_y);

    for (auto & elm : Y) cout << elm << ' ';
    cout << endl;

    // method.compute_truncated(Y, )


    // frlap::gdm::far_field_estimator<double, far_field::zero> ffestim_0();
    // frlap::gdm::far_field_estimator<double, far_field::general> ffestim_g(
    //     y, ealpha, deltax, a, b);
    // frlap::gdm::far_field_estimator<double, far_field::algebraic> ffestim_a(
    //     y, ealpha, deltax, a, b, a0, b0, edecay);

    if (true) {
        cout << "Writing...\n";
    }

    return EXIT_SUCCESS;
}

// Write the results to a file:
void write_results(std::string filename,
                   std::size_t n0,
                   std::size_t ja,
                   double const * const X,
                   double const * const Y,
                   double const * const FLY_0,
                   double const * const FF,
                   double const * const Exact) {
    std::FILE *fp;

    if ((fp = std::fopen(filename.c_str(), "w")) == NULL) {
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
                            FLY_0[j]+FF[j],
                            Exact[j+ja],
                            fabs(Exact[j+ja]- (FLY_0[j] + FF[j])));
    std::fclose(fp);
}
