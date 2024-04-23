// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <cstddef>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include "./simple_calculator.hpp"

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

int main(int argc, char * argv[]) {
    using std::cout;
    using std::endl;
    
    // Input parameters
    CoefficientsKind coefficients_kind;
    unsigned         domain_multiplier = 3;
    unsigned         grid_refinement   = 1;


    (void) std::puts("\n\t*** Benchmark 7.0.1 ***\n");

    // Assign input
    if (argc > 1) {
        try {
            coefficients_kind = CoeffKindFromStr.at(argv[1]);
            cout << CoeffKindMessage.at(coefficients_kind)
                 << endl;
        } catch (std::out_of_range const & e) {
            cout << "invalid input parameter. ";
            std::string const what{e.what()};
            if (what == "map::at") {
                coefficients_kind = CoefficientsKind::HuangOberman2;
                cout << CoeffKindMessage.at(coefficients_kind)
                     << endl;
            } else {
                cout << "error: " << e.what() << endl; 
            }
        }
    }
    if (argc > 2) {
        domain_multiplier = std::atoi(argv[2]);
    }
    if (argc > 3) {
        grid_refinement   = std::atoi(argv[3]);
    }

    // TODO(gffrnl): Verificação de erros de params. de entrada


    // Computes the execution time
    auto time_0 = std::chrono::high_resolution_clock::now();
    
    double const ealpha = 0.4;
    //double const deltax = 0.1;
    double const deltax = 0.08;
    //double const a = -6.0;
    double const a = -12.0;
    //double const b = +6.0;
    double const b = +12.0;
    //double const a0 = -2.0;
    double const a0 = -4.0;
    //double const b0 = +2.0;
    double const b0 = +4.0;
    //double const n = 121;
    double const n = 31;
    double const edecay = 1-ealpha;

    double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
    
    auto y       = [&ealpha, &sgn](double const & x) -> double {
        return sgn * pow(1.0 + x*x, -0.5 + ealpha/2.0);
    };
    auto frLap_y = [&ealpha, &sgn](double const & x) -> double {
        return sgn * pow(1.0 + x*x, -0.5 - ealpha/2.0) * exp2(ealpha)
            * tgamma(0.5+ealpha/2.0) / tgamma(0.5-ealpha/2.0);
    };

    b118::grid<double> G({a, b}, n);
    b118::grid<double> G0 = G.subgrid({a0, b0});

    b118::grid_function<double> Y   (G, y);
    b118::grid_function<double> Y0  (G0, y);
    b118::grid_function<double> FLYe(G0, frLap_y);
    b118::grid_function<double> FLY0(G0);
    b118::grid_function<double> FF  (G0);
    b118::grid_function<double> FLY (G0);

    // frlap::generalized_differences<double> method(ealpha, G);
    frlap::generalized_differences<
      double,
      b118::frlap::gdm::coefficients::huang_oberman_2
      > method(ealpha, G);

    auto time_1 = std::chrono::high_resolution_clock::now();
    auto delta_1 =
            std::chrono::duration_cast<std::chrono::milliseconds>(time_1 - time_0);
    cout << "Duration of the initialization: "
         << delta_1.count() << "ms" << endl;

    auto time_2 = std::chrono::high_resolution_clock::now();
    method.compute_truncated(Y, &FLY0);
    auto time_3 = std::chrono::high_resolution_clock::now();
    auto delta_2 =
            std::chrono::duration_cast<std::chrono::milliseconds>(time_3 - time_2);
      cout << "Duration of 'compute_truncated': " << delta_2.count() << "ms" << endl;


    auto time_4 = std::chrono::high_resolution_clock::now();
    // method.compute_far_field(far_field::general(y), &FF);
    // method.compute_far_field(far_field::algebraic(edecay, y(a0), y(b0)), &FF); 
    // method.compute_far_field(far_field::algebraic(edecay,
    //                                               y(G0.lendpoint()),
    //                                               y(G0.rendpoint())), &FF); 
    auto time_5 = std::chrono::high_resolution_clock::now();
    auto delta_3 =
            std::chrono::duration_cast<std::chrono::milliseconds>(time_5 - time_4);
      cout << "Duration of 'compute_far_field': " << delta_3.count() << "ms" << endl;
    
    

    if (true) {
        cout << "Writing...\n";
        auto time_6 = std::chrono::high_resolution_clock::now();
        write_results("bench-701.dat", G0, Y0, FLY0, FF, FLYe);
        auto time_7 = std::chrono::high_resolution_clock::now();
        auto delta_4 =
            std::chrono::duration_cast<std::chrono::milliseconds>(time_7 - time_6);
      cout << "Duration to write results to the file: " << delta_4.count() << "ms" << endl;
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
                            std::fabs(FLYe[j] - value));
    }
    std::fclose(fp);
}
