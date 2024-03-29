// Copyright 2014 Guilherme F. Fornel <gffrnl@gmail.com>
// Compilar com:
// g++ -std=c++17 main.cpp -Wall -Wpedantic -lfftw3f -lfftw3 -lfftw3l -lfftw3q -lopenblas -DNDEBUG

#include <iostream>
#include <chrono>  // NOLINT
#include <b118/linspace.hpp>
#include <b118/closest_sorted.hpp>
#include <b118/fftw/managed.hpp>
#include <b118/frlap/gdm/coefficients/huang_oberman_1.hpp>
#include <b118/frlap/gdm/coefficients/huang_oberman_2.hpp>
#include "./generalized_differences.hpp"

using real_t = double;

namespace Benchmark701 {
  real_t      ealpha =  0.4;  // fractional Laplacian exponent
  real_t      a0     = -2.0;  // left inner endpoint
  real_t      b0     = +2.0;  // right inner endpoint

  real_t y(real_t x) {
      real_t const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
      return sgn * std::pow(1.0 + x*x, -0.5 + ealpha/2.0);
  }

  real_t frLap_y(real_t x) {
      real_t const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
      return sgn * std::exp2(ealpha) *
          std::tgamma(0.5+ealpha/2.0) / std::tgamma(0.5-ealpha/2.0) *
          std::pow(1.0 + x*x, -0.5-ealpha/2.0);
  }

  void compute_far_field(std::vector<real_t> const & x ,
                         std::size_t                 ja,
                         std::size_t                 jb,
                         std::vector<real_t>       & ff) {  // NOLINT
      using namespace b118::frlap::gdm;  // NOLINT
      namespace far_field = b118::frlap::gdm::far_field;
      real_t const deltax = x.at(2) - x.at(1);
      far_field_estimator<real_t, far_field::general> ffest(y, ealpha, deltax,
        x.front(), x.back());
      for (std::size_t j = ja; j <= jb; ++j) {
        ff[j-ja] = ffest(x[j]);
      }
  }
}

int main(int argc, char * argv[]) {
    using std::cout;
    using std::clog;
    using std::cerr;
    using std::endl;
    using b118::exclude;
    using Benchmark701::ealpha;
    using Benchmark701::a0;
    using Benchmark701::b0;

    namespace frlap = b118::frlap;
    namespace coefficients = b118::frlap::gdm::coefficients;
    namespace far_field = b118::frlap::gdm::far_field;

  // params
  std::string tipo_coefficientes = argc > 1? argv[1] : "huob2";
  int multiplicador_de_dominio = argc > 2? std::atoi(argv[2]) : 3;
  int refino_da_malha = argc > 3? std::atoi(argv[3]) : 1;

  if (!(tipo_coefficientes == "spec")      &&
      !(tipo_coefficientes == "spec_qawo") &&
      !(tipo_coefficientes == "spec_thsh") &&
      !(tipo_coefficientes == "gormai")    &&
      !(tipo_coefficientes == "huob1")     &&
      !(tipo_coefficientes == "huob2")     &&
      !(tipo_coefficientes == "c3point")) {
    cerr << "primeiro argumento deve ser:\n"
            "  `spec_qawo`: coeficientes espectrais\n"
            "  `spec_qawo`: coeficientes espectrais (qawo)\n"
            "  `spec_thsh`: coeficientes espectrais (tanh_sinh)\n"
            "  `gormai`   : coeficientes de Gorenflo & Mainardi\n"
            "  `huob1`    : coeficientes de Huang & Oberman linear\n"
            "  `huob2`    : coeficientes de Huang & Oberman quadratico\n"
            "  `c3point`  : coeficientes da periodizacao da regra de 3 pontos"
         << endl;
    return EXIT_FAILURE;
  }

  real_t deltax;
  std::vector<real_t> unigrid;  // Storage for the uniform grid
  std::vector<real_t> fnyvals;  // Storage for funciton y values

  real_t a = a0 * multiplicador_de_dominio;
  real_t b = b0 * multiplicador_de_dominio;
  std::size_t n = 40 * multiplicador_de_dominio * refino_da_malha + 1;

  (void) std::puts("\n\t*** Benchmark 7.0.1 ***\n");

  if (tipo_coefficientes == "spec") {
    cout << "[using spectral]" << endl;
  } else if (tipo_coefficientes == "spec_qawo") {
    cout << "[using spectral_qawo]" << endl;
  } else if (tipo_coefficientes == "spec_thsh") {
    cout << "[using spectral_thsh]" << endl;
  } else if (tipo_coefficientes == "gormai") {
    cout << "[using gorenflo_mainardi]" << endl;
  } else if (tipo_coefficientes == "huob1") {
    cout << "[using huang_oberman_linear]" << endl;
  } else if (tipo_coefficientes == "huob2") {
    cout << "[using huang_oberman_quadratico]" << endl;
  } else if (tipo_coefficientes == "c3point") {
    cout << "[using centered_3_point_periodized]" << endl;
  } else {
    cerr << "primeiro argumento invalido. abortando..." << endl;
    abort();
  }
  cout << "[using multiplicador_de_dominio = " << multiplicador_de_dominio
        << ']' << endl;
  cout << "[using refino_da_malha = " << refino_da_malha << ']' << endl;
  cout << "[using a = " << a << ']' << endl;
  cout << "[using b = " << b << ']' << endl;
  cout << "[using n = " << n << ']' << endl;

  // Calcula o tempo de execução
  auto tempo_0 = std::chrono::high_resolution_clock::now();

  unigrid = b118::linspace<std::vector, real_t>(a, b, n,
                                                exclude::none,
                                                &deltax);
  cout << "\n deltax = " << deltax << endl;
  assert(unigrid.size() == n);

#ifndef NDEBUG
  cout << "uniform grid: {" << n << ';' << deltax << "} (";
  for (auto & x : unigrid)
    cout << x << ' ';
  cout << "\b)" << endl;
#endif

  fnyvals.resize(n);
  for (std::size_t k{0}; k < n; ++k) {
    fnyvals.at(k) = Benchmark701::y(unigrid.at(k));
  }

#ifndef NDEBUG
  cout << "values of function y at grid points: (";
  for (auto & yk : fnyvals)
    cout << yk << ' ';
  cout << "\b)" << endl;
#endif

  std::size_t ja = b118::closest_sorted(unigrid.cbegin(), unigrid.cend(),
                                        a0);
  std::size_t jb = b118::closest_sorted(unigrid.cbegin(), unigrid.cend(),
                                        b0);

  assert(ja < jb);

#ifndef NDEBUG
  cout << "ja: " << ja << endl;
  cout << "jb: " << jb << endl;
#endif
  auto tempo_1 = std::chrono::high_resolution_clock::now();
  auto delta_1 =
    std::chrono::duration_cast<std::chrono::milliseconds>(tempo_1 - tempo_0);
  cout << "Duracao da inicializacao: " << delta_1.count() << "ms" << endl;

  std::vector<real_t> FLY0;
  decltype(tempo_1) tempo_2;
  decltype(tempo_1) tempo_3;

  if (tipo_coefficientes == "spec") {
    std::abort();
    // b118::frlap::gdm::trunc_uniform_spec method(ealpha, deltax);
    // tempo_2 = std::chrono::high_resolution_clock::now();
    // method.compute(fnyvals, ja, jb, FLY0);
    // tempo_3 = std::chrono::high_resolution_clock::now();
  } else if (tipo_coefficientes == "spec_qawo") {
    std::abort();
    // b118::frlap::gdm::trunc_uniform_spec_qawo method(ealpha, deltax);
    // tempo_2 = std::chrono::high_resolution_clock::now();
    // method.compute(fnyvals, ja, jb, FLY0);
    // tempo_3 = std::chrono::high_resolution_clock::now();
  } else if (tipo_coefficientes == "spec_thsh") {
    std::abort();
    // b118::frlap::gdm::trunc_uniform_spec_thsh method(ealpha, deltax);
    // tempo_2 = std::chrono::high_resolution_clock::now();
    // method.compute(fnyvals, ja, jb, FLY0);
    // tempo_3 = std::chrono::high_resolution_clock::now();
  } else if (tipo_coefficientes == "gormai") {
    std::abort();
    // b118::frlap::gdm::trunc_uniform_gormai method(ealpha, deltax);
    // tempo_2 = std::chrono::high_resolution_clock::now();
    // method.compute(fnyvals, ja, jb, FLY0);
    // tempo_3 = std::chrono::high_resolution_clock::now();
  } else if (tipo_coefficientes == "huob1") {
    frlap::generalized_differences<
            real_t,
            coefficients::huang_oberman_1,
            far_field::general
            > method(ealpha, a, b, n);
    tempo_2 = std::chrono::high_resolution_clock::now();
    method.compute_0(fnyvals, ja, jb, FLY0);
    tempo_3 = std::chrono::high_resolution_clock::now();
  } else if (tipo_coefficientes == "huob2") {
    frlap::generalized_differences<
            real_t,
            coefficients::huang_oberman_2,
            far_field::general
            > method(ealpha, a, b, n);
    tempo_2 = std::chrono::high_resolution_clock::now();
    method.compute_0(fnyvals, ja, jb, FLY0);
    tempo_3 = std::chrono::high_resolution_clock::now();
  } else if (tipo_coefficientes == "c3point") {
    std::abort();
    // cout << "[using centered_3_point_periodized]" << endl;
    // b118::frlap::gdm::trunc_uniform_c3point method(ealpha, deltax);
    // tempo_2 = std::chrono::high_resolution_clock::now();
    // method.compute(fnyvals, ja, jb, FLY0);
    // tempo_3 = std::chrono::high_resolution_clock::now();
  } else {
    cerr << "primeiro argumento invalido. abortando..." << endl;
    std::abort();
  }
  assert(FLY0.size() == (jb - ja + 1));

  auto delta_2 =
    std::chrono::duration_cast<std::chrono::milliseconds>(tempo_3 - tempo_2);
  cout << "Duração do 'compute': " << delta_2.count() << "ms" << endl;

#ifndef NDEBUG
  cout << "FLY0: (";
  for (auto & elem : FLY0)
      cout << elem << ' ';
  cout << "\b)" << endl;
#endif

  std::size_t n0 = FLY0.size();  // == jb - ja + 1

  std::vector<real_t> FF(n0);
  Benchmark701::compute_far_field(unigrid, ja, jb, FF); // NOLINT
  auto tempo_4 = std::chrono::high_resolution_clock::now();
  auto delta_3 =
    std::chrono::duration_cast<std::chrono::milliseconds>(tempo_4 - tempo_3);
  cout << "Duração do 'compute_far_field': " << delta_3.count() << "ms" << endl;

  // Calculate error
  real_t sum_error = 0;
  real_t max_error = 0;
  real_t sum_abs_error = 0;
  for (std::size_t j = 0; j < n0; j++) {
    real_t error = Benchmark701::frLap_y(unigrid[j+ja]) - (FLY0[j]+FF[j]);
    sum_error += error;
    max_error = std::max(max_error, std::fabs(error));
    sum_abs_error += std::fabs(error);
  }
  std::cout << "\n------\n\n";
  std::cout << std::scientific;
  // std::cout << "average error = " << sum_error/n0 << '\n';
  std::cout << "average abs error = " << sum_abs_error/n0 << '\n';
  std::cout << "max error = " << max_error << '\n';

  // Write the results to a file:
  if (true) {
    std::FILE *fp;

    if ((fp = std::fopen("bench-701.dat", "w")) == NULL) {
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
                          unigrid[j+ja], fnyvals[j+ja], FLY0[j], FF[j],
                          FLY0[j]+FF[j],
                          Benchmark701::frLap_y(unigrid[j+ja]),
                          std::fabs(Benchmark701::frLap_y(unigrid[j+ja])-
                                    (FLY0[j] + FF[j])));
    std::fclose(fp);
  }
    return EXIT_SUCCESS;
}
