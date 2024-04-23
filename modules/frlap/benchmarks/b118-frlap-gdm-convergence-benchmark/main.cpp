// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#include <cstdlib>
#include <memory>
#include <string>
#include <map>
#include <chrono>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <b118/grid.hpp>
#include <b118/grid_function.hpp>
#include <b118/frlap/generalized_differences.hpp>
#include <b118/frlap/gdm/coefficients.hpp>
#include <b118/frlap/gdm/far_field.hpp>
#include "./ConvergenceBenchmark.hpp"

int const nmin = 2;
int const mmin = 0;
char const * input_err_msg = "*** input error: ";

static std::map<int, int> L;
static std::map<int, double> h;

std::string label_from_input(double alpha, int nmax, int mmax);


template<typename Real>
Real foo(ConvergenceBenchmark::Problem_ptr<double> prob_ptr,
         int n, int m,
         std::chrono::milliseconds durs[2]);

template<typename Real>
Real bar(ConvergenceBenchmark::Problem_ptr<Real> prob_ptr,
         b118::grid<Real> GhJ, b118::grid<Real> GhJ0,
         std::chrono::milliseconds durs[2]);

template<
    typename Real,
    template<typename> class CoeffGenerator
    >
Real baz(ConvergenceBenchmark::Problem_ptr<Real> prob_ptr,
         b118::grid<Real> GhJ, b118::grid<Real> GhJ0,
         std::chrono::milliseconds durs[2]);


int main(int argc, char * argv[]) {
    using std::cout;
    using std::cerr;
    using std::endl;
    namespace method = b118::frlap::gdm::coefficients;
    namespace far_field = b118::frlap::gdm::far_field;
    
    int nmax;
    int mmax;

    // Pass input parameters
    if (argc <= 2) {
        cout << " Usage:\n Example nmax mmax" << endl;
        return EXIT_FAILURE;
    } else {
        nmax  = std::atoi(argv[1]);
        mmax  = std::atoi(argv[2]);
    }


    ConvergenceBenchmark::Problem_ptr<double> prob_ptr
            = ConvergenceBenchmark::create_problem<double>();

    // TODO(gffrnl): check alpha separately;

    
    // Check input parameters 
    {
        bool has_input_error = false;
        if (prob_ptr->alpha < 0 || prob_ptr->alpha > 2) {
            cerr << input_err_msg
                 << "alpha must lie in the range [0, 2]"
                 << endl;
            has_input_error = true;
        }
        if (nmax < nmin) {
            cerr << input_err_msg
                 << "nmax must be greater than or equal to 3"
                 << endl;
            has_input_error = true;
        }
        if (mmax < mmin) {
            cerr << input_err_msg
                 << "nmax must be greater than or equal to 0"
                 << endl;
            has_input_error = true;
        }
        if (has_input_error) {
            cerr << "exiting ..." << endl;
            return EXIT_FAILURE;
        }
    }

    // Print input parameters
    cout << "alpha = " << prob_ptr->alpha << '\n';
    cout << "nmax  = " << nmax  << '\n';
    cout << "mmax  = " << mmax  << '\n';

    std::string const label = label_from_input(prob_ptr->alpha, nmax, mmax);
    cout << "label = " << label << '\n';

    for (int n = nmin; n <= nmax; ++n)
        L[n] = 1 << n;
    for (int m = mmin; m <= mmax; ++m)
        h[m] = 4.0 / (5 * (1 << m));

    std::string const outfname{prob_ptr->get_label() + "_" + label +
                               ".dat"};
    std::string const timifname{prob_ptr->get_label() + "_" + label +
                               "_tim_init.dat"};
    std::string const timcfname{prob_ptr->get_label() + "_" + label +
                               "_tim_comp.dat"};
    
    std::fstream outfile;
    outfile.open(outfname, std::ios::out);
    if (outfile.is_open()) {
        outfile.precision(8);
        outfile << std::setw(14) << std::right << "h";
        for (int n = nmin; n <= nmax; ++n) {
            outfile << '\t';
            outfile << std::setw(14) << std::right
                    << (std::string("L=") + std::to_string(L[n]));
        }
        outfile << endl;
    } else {
        std::clog << "!!! Could not open file " << outfname << "to write !!!"
                  << std::endl;
    }

    std::fstream timifile;
    timifile.open(timifname, std::ios::out);
    if (timifile.is_open()) {
        timifile.precision(8);
        timifile << std::setw(14) << std::right << "h";
        for (int n = nmin; n <= nmax; ++n) {
            timifile << '\t';
            timifile << std::setw(14) << std::right
                    << (std::string("L=") + std::to_string(L[n]));
        }
        timifile << endl;
    } else {
        std::clog << "!!! Could not open file " << timifname
                  << " to write !!!"
                  << std::endl;
    }

    std::fstream timcfile;
    timcfile.open(timcfname, std::ios::out);
    if (timcfile.is_open()) {
        timcfile.precision(8);
        timcfile << std::setw(14) << std::right << "h";
        for (int n = nmin; n <= nmax; ++n) {
            timcfile << '\t';
            timcfile << std::setw(14) << std::right
                    << (std::string("L=") + std::to_string(L[n]));
        }
        timcfile << endl;
    } else {
        std::clog << "!!! Could not open file " << timcfname
                  << " to write !!!"
                  << std::endl;
    }

    for (int m = mmin; m <= mmax; ++m) {
        if (outfile.is_open()) {
            outfile << std::scientific << h[m] << " ";
        } else {
            std::clog << "!!! output file " << outfname << " not opened !!!"
                      << std::endl;
        }
        if (timifile.is_open()) {
            timifile << std::scientific << h[m] << " ";
        } else {
            std::clog << "!!! times file " << timifname << " not opened !!!"
                      << std::endl;
        }
        if (timcfile.is_open()) {
            timcfile << std::scientific << h[m] << " ";
        } else {
            std::clog << "!!! times file " << timcfname << " not opened !!!"
                      << std::endl;
        }

        for (int n = nmin; n <= nmax; ++n) {
            double max_abs_err;
            std::chrono::milliseconds durs[2];

            max_abs_err = foo<double>(prob_ptr, n, m, durs);

            cout << "Duration of the initialization: "
                 << durs[0].count() << "ms" << endl;
            cout << "Duration of 'compute': "
                 << durs[1].count() << "ms" << endl;
            cout << "||flye - fly ||max  = "
                 << max_abs_err << endl;
            
            if (outfile.is_open()) {
                outfile << '\t';
                outfile << std::scientific << max_abs_err << " ";
            } else {
                std::clog << "!!! output file " << outfname
                          << " not opened !!!"
                          << std::endl;
            }
            if (timifile.is_open()) {
                timifile << '\t';
                timifile << std::setw(14) << std::right << durs[0].count();
                timifile << " ";
            } else {
                std::clog << "!!! times file " << timifname
                          << " not opened !!!"
                          << std::endl;
            }
            if (timcfile.is_open()) {
                timcfile << '\t';
                timcfile << std::setw(14) << std::right << durs[1].count();
                timcfile << " ";
            } else {
                std::clog << "!!! times file " << timcfname
                          << " not opened !!!"
                          << std::endl;
            }
        }

        if (outfile.is_open()) {
            outfile << '\n';
        } else {
            std::clog << "!!! output file " << outfname << " not opened !!!"
                      << std::endl;
        }
        if (timifile.is_open()) {
            timifile << '\n';
        } else {
            std::clog << "!!! tim file " << timifname << " not opened !!!"
                      << std::endl;
        }
        if (timcfile.is_open()) {
            timcfile << '\n';
        } else {
            std::clog << "!!! tim file " << timcfname << " not opened !!!"
                      << std::endl;
        }
    }

    if (outfile.is_open() ) outfile.close();
    if (timifile.is_open()) timifile.close();
    if (timcfile.is_open()) timcfile.close();

    std::cout << std::endl;
    std::cout << "OUTPUT    FILE SAVED AS " << outfname << std::endl;
    std::cout << "DURATIONS FILE SAVED AS " << timifname << std::endl;
    std::cout << "DURATIONS FILE SAVED AS " << timcfname << std::endl;
    
    return EXIT_SUCCESS;
}


std::string label_from_input(double alpha, int nmax, int mmax) {
    std::string alph_str("a");
    std::string nmax_str("n");
    std::string mmax_str("m");
    {
        std::ostringstream oss;
        oss.precision(4);
        oss << std::fixed << alpha;
        alph_str += oss.str();
        auto it = std::remove(alph_str.begin(), alph_str.end(), '.');
        alph_str = std::string(alph_str.begin(), it);
    }
    nmax_str += std::to_string(nmax);
    mmax_str += std::to_string(mmax);
    return alph_str + nmax_str + mmax_str;
}



template<typename Real>
Real foo(ConvergenceBenchmark::Problem_ptr<double> prob_ptr, int n, int m,
         std::chrono::milliseconds durs[2]) {
    std::size_t const J  = 15 * (1 << (n+m-1)) + 1;
    std::size_t const J0 =  5 * (1 << (n+m-1)) + 1;

    b118::grid<Real> const GhJ ({-3.0*L[n], 3.0*L[n]}, J );
    b118::grid<Real> const GhJ0({-1.0*L[n], 1.0*L[n]}, J0);

#ifndef NDEBUG
    if (!(GhJ0.is_subgrid(GhJ))) {
        throw std::runtime_error("!(GhJ0.is_subgrid(GhJ))");
    }
#endif

    return bar(prob_ptr, GhJ, GhJ0, durs);
}


template<typename Real>
Real bar(ConvergenceBenchmark::Problem_ptr<Real> prob_ptr,
         b118::grid<Real> GhJ, b118::grid<Real> GhJ0,
         std::chrono::milliseconds durs[2])
{
    using namespace b118::frlap::gdm::coefficients;
    
    Real max_abs_err = static_cast<Real>(0);
    switch (prob_ptr->method) {
        case (ConvergenceBenchmark::Method::spectral):
            max_abs_err = baz<Real, spectral>(prob_ptr, GhJ, GhJ0, durs);
            break;

        case (ConvergenceBenchmark::Method::spectral_qawo):
            max_abs_err = baz<Real, spectral_qawo>(prob_ptr, GhJ, GhJ0, durs);
            break;

        case (ConvergenceBenchmark::Method::gorenflo_mainardi):
            max_abs_err = baz<Real, gorenflo_mainardi>(prob_ptr, GhJ, GhJ0, durs);
            break;

        case (ConvergenceBenchmark::Method::huang_oberman_1):
            max_abs_err = baz<Real, huang_oberman_1>(prob_ptr, GhJ, GhJ0, durs);
            break;

        case (ConvergenceBenchmark::Method::huang_oberman_2):
            max_abs_err = baz<Real, huang_oberman_2>(prob_ptr, GhJ, GhJ0, durs);
            std::cout << "Using coefficients huang_oberman_2" << std::endl;
            break;

        default:  // ConvergenceBenchmark::Method::cper3point
            max_abs_err = baz<Real, cper3point>(prob_ptr, GhJ, GhJ0, durs);
    }

    return max_abs_err;

}

template<
    typename Real,
    template<typename> class CoeffGenerator
    >
Real baz(ConvergenceBenchmark::Problem_ptr<Real> prob_ptr,
         b118::grid<Real> GhJ, b118::grid<Real> GhJ0,
         std::chrono::milliseconds durs[2])
{
    using namespace std::chrono;
    using clock = high_resolution_clock;
    using ms = milliseconds;
    namespace far_field = b118::frlap::gdm::far_field;
    namespace FarField = ConvergenceBenchmark::FarField;

    time_point<clock> time[4];

    time[0] = clock::now();

    Real const alpha = prob_ptr->alpha;
    auto y       = [&prob_ptr](Real const & x) -> Real {
        return prob_ptr->y(x);
    };
    auto frLap_y = [&prob_ptr](Real const & x) -> Real {
        return prob_ptr->frLap_y(x);
    };

    b118::frlap::generalized_differences<Real, CoeffGenerator> gdm(alpha, GhJ);
    b118::grid_function<Real> fly(GhJ0);
        
    if (std::holds_alternative<FarField::algebraic<Real>>(prob_ptr->farfield)) {
        std::cout << "using farfield algebraic decay" << std::endl;
        std::cout << "prob_ptr->get_decay() = "
                  <<  prob_ptr->get_decay()
                  << std::endl;
        time[1] = clock::now();
        gdm.compute(far_field::algebraic(prob_ptr->get_decay(),
                                         y(GhJ0.lendpoint()),
                                         y(GhJ0.rendpoint())),
                    b118::grid_function<Real>{GhJ, y}, &fly);
        time[2] = clock::now();
    } else {
        switch (std::get<unsigned>(prob_ptr->farfield)) {
            case FarField::zero:
                time[1] = clock::now();
                gdm.compute_truncated(b118::grid_function<Real>{GhJ, y}, &fly);
                time[2] = clock::now();
                break;
            default:  // FarField::general
                time[1] = clock::now();
                gdm.compute(far_field::general(y),
                            b118::grid_function<Real>{GhJ, y}, &fly);
                time[2] = clock::now();
                break;
        }
    }

    durs[0] = duration_cast<ms>(time[1] - time[0]);
    durs[1] = duration_cast<ms>(time[2] - time[1]);

    return b118::distance::max(b118::grid_function<Real>{GhJ0, frLap_y},
                               fly);
}
