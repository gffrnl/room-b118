// Compilar com g++ -std=c++17 cspectral.cpp -o cspectral -lgsl
#include <gsl/gsl_integration.h>
#include <cstddef>
#include <cmath>
#include <cassert>
#include <iostream>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::constants;


static double f(double x, void* p) {
  double const expon = *static_cast<double *>(p);
  return std::pow(x, expon);
}

gsl_integration_workspace  * w;
gsl_integration_qawo_table * t;
gsl_function F;
std::size_t levels = 30;


double com_qawo(std::size_t k, double * err = nullptr) {
  double err1;
  double ret;
  gsl_integration_qawo_table_set(t, k, M_PI, GSL_INTEG_COSINE);
  gsl_integration_qawo(&F, 0.0, 1.1e-12, 1.1e-12, levels, w, t,
                       &ret,
                       &err1);
  if (err != nullptr) *err = err1;
  return ret;
}

template<typename Real>
Real serie_ascendente(Real ealpha, std::size_t k) {
  assert(k > 0);
  Real x = static_cast<Real>(k*k)*pi_sqr<Real>()/static_cast<Real>(4);
  Real x2 = x*x;
  Real C1 = static_cast<Real>(2) - static_cast<Real>(4)/(ealpha + static_cast<Real>(3)) * x;
  Real C2 = C1 * (static_cast<Real>(2)/static_cast<Real>(3)
                  - static_cast<Real>(4)/(static_cast<Real>(3)*ealpha+static_cast<Real>(15)));
  for (std::size_t n = 1; n < 30; n += 2) {
    C1 += C1 * (
            (static_cast<Real>(1) - static_cast<Real>(4)/(ealpha + static_cast<Real>(2*n) + static_cast<Real>(5))
            +
             static_cast<Real>(1) / (static_cast<Real>(n*n)+static_cast<Real>(2*n)+static_cast<Real>(3)/static_cast<Real>(4))
             ) * x2
    );
  }
  for (std::size_t n = 2; n < 30; n += 2) {
    C2 += C2 * (
            (static_cast<Real>(1) - static_cast<Real>(4)/(ealpha + static_cast<Real>(2*n) + static_cast<Real>(5))
            +
             static_cast<Real>(1) / (static_cast<Real>(n*n)+static_cast<Real>(2*n)+static_cast<Real>(3)/static_cast<Real>(4))
             ) * x2
    );
  }

  
  return C2-C1;
}

int main() {
  using std::cout;
  using std::endl;
  double const ealpha = 0.4;
  double const cpiealpha = (ealpha + 1) / std::pow(M_PI, ealpha + 1);
  
  w = gsl_integration_workspace_alloc(levels);
  t = gsl_integration_qawo_table_alloc(0, M_PI, GSL_INTEG_COSINE, levels);
  F.function = &f;
  F.params = static_cast<void *>(const_cast<double *>(&ealpha));

  for (std::size_t i = 1; i < 30; i += 10) {
    double qawo_est_err = 99;
    cout << i << '\t'
         << cpiealpha * com_qawo(i, &qawo_est_err) << "\t[" << qawo_est_err << "]\t\t"
         << serie_ascendente<double>(ealpha, i) << '\t'
         << endl;
  }
  
  gsl_integration_qawo_table_free(t);
  gsl_integration_workspace_free(w);  
  return 0;
}
