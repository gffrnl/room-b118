// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <cstddef>
#include <cmath>
#include <iostream>
#include <string>
#include <b118/linspace.hpp>
#include <b118/frlap/generalized_differences.hpp>
#include <b118/frlap/generalized_differences_methods/coefficients/huang_oberman_2.hpp>
#include <b118/frlap/generalized_differences_methods/far_field.hpp>


#include <b118/grid.hpp>

void write_results(std::string                 filename,
                   b118::grid<double>          const & G0,
                   b118::grid_function<double> const & Y0,
                   b118::grid_function<double> const & FLY0,
                   b118::grid_function<double> const & FF,
                   b118::grid_function<double> const & FLYe);

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

    b118::grid G({a, b}, n);
    b118::grid G0 = G.subgrid({a0, b0});

    b118::grid_function Y   (G, y);
    b118::grid_function Y0  (G0, y);
    b118::grid_function FLYe(G0, frLap_y);
    b118::grid_function FLY0(G0);
    b118::grid_function FF  (G0);

    // cout << "Y = ";
    // for (std::size_t k = 0; k < G.numnodes(); ++k)
    //     cout << Y[k] << ", ";
    // cout << endl;
    // cout << "FLYe = ";
    // for (std::size_t k = 0; k < G0.numnodes(); ++k)
    //     cout << FLYe[k] << ", ";
    // cout << endl;
    // cout << "FLY0 = ";
    // for (std::size_t k = 0; k < G0.numnodes(); ++k)
    //     cout << FLY0[k] << ", ";
    // cout << endl;

    
    frlap::generalized_differences method(ealpha, {a, b}, n);

    method.compute_truncated(Y, &FLY0); // TODO(gffrnl): error if the grid of FLY0 is not subgrid of the grid of Y
    // cout << "FLY0 = ";
    // for (std::size_t k = 0; k < G0.numnodes(); ++k)
    //     cout << FLY0[k] << ", ";
    // cout << endl;

    // // Zero far-fieldestimator:
    // frlap::gdm::far_field_estimator<double, far_field::zero>
    // ffestim;
    
    // General far-field estimator:
    frlap::gdm::far_field_estimator<double, far_field::general>
    ffestim(y, method.ealpha(), G);
    // TODO(gffrnl): Trocar para ffestim(y, method.ealpha(), method.get_grid()) ??;
    
    // // Algebric decay far-field estimator:
    // frlap::gdm::far_field_estimator<double, far_field::algebraic>
    // ffestim(y, method.ealpha(), G, G0, edecay);
     
    // for (std::size_t k = 0; k < G0.numnodes(); ++k)
    //     FF[k] = ffestim(G0[k]);

    method.compute_far_field(ffestim, &FF);
    
    // cout << "FF = ";
    // for (std::size_t k = 0; k < G0.numnodes(); ++k)
    //     cout << FF[k] << ", ";
    // cout << endl;

    
    
    // frlap::gdm::far_field_estimator<double, far_field::zero> ffestim_0();
    // frlap::gdm::far_field_estimator<double, far_field::general> ffestim_g(
    //     y, ealpha, deltax, a, b);
    // frlap::gdm::far_field_estimator<double, far_field::algebraic> ffestim_a(
    //     y, ealpha, deltax, a, b, a0, b0, edecay);

    if (true) {
        cout << "Writing...\n";
        write_results("bench-701.dat", G0, Y0, FLY0, FF, FLYe);
    }

    return EXIT_SUCCESS;
}

// Write the results to a file:
void write_results(std::string                 filename,
                   b118::grid<double>          const & G0,
                   b118::grid_function<double> const & Y0,
                   b118::grid_function<double> const & FLY0,
                   b118::grid_function<double> const & FF,
                   b118::grid_function<double> const & FLYe) {
    std::FILE *fp;

    if ((fp = std::fopen(filename.c_str(), "w")) == NULL) {
        (void) std::fprintf(stderr,
                            "error: fopen() :"
                            "cannot open file to write\n");
        std::abort();
    }

    (void) std::fprintf(fp,
                        "%-22s\t%-22s\t%-22s\t%-22s\t%-22s\t%-22s\t%-22s\n",
                        "X", "Y", "FLY0", "FF", "FLY", "Exact", "Error");
    for (std::size_t j = 0; j < G0.numnodes(); ++j) {
        auto const value = FLY0[j] + FF[j];
        (void) std::fprintf(fp,
                            "%+22.15E\t%+22.15E\t%+22.15E\t%+22.15E\t"
                            "%+22.15E\t%+22.15E\t%+22.15E\n",
                            G0[j], Y0[j], FLY0[j], FF[j],
                            FLY0[j]+FF[j],
                            FLYe[j],
                            fabs(FLYe[j]- value));
    }
    std::fclose(fp);
}
