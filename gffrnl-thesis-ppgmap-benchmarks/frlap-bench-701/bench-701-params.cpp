// Copyright 2024 Guilherme F. Fornel

#include <cstddef>
#include <cmath>
#include <cassert>
#include <vector>
#include <b118/frlap.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>


namespace Benchmark701 {
  double      ealpha =  0.4;  // fractional Laplacian exponent
  double      a0     = -2.0;  // left inner endpoint
  double      b0     = +2.0;  // right inner endpoint

  double y(double x) {
      double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
      return sgn * std::pow(1.0 + x*x, -0.5 + ealpha/2.0);
  }

  double frLap_y(double x) {
      extern double ealpha;
      double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
      return sgn * std::exp2(ealpha) *
          std::tgamma(0.5+ealpha/2.0) / std::tgamma(0.5-ealpha/2.0) *
          std::pow(1.0 + x*x, -0.5-ealpha/2.0);
  }


  // void compute_far_field(std::vector<double> const & x ,
  //                        std::size_t                 ja,
  //                        std::size_t                 jb,
  //                        std::vector<double>       & ff)  { // NOLINT
  //     double C1a = b118::frlap::normal_const<1>(ealpha);
  //     double ebeta = 1.0 - ealpha;
  //     double h = x[1]-x[0];

  //     // for (std::size_t j = ja; j <= jb; ++j) {
      //     ff[j-ja] = -C1a/(ealpha+ebeta) * (
      //     y(x[ja])*std::pow(std::fabs(x[ja]),beta)/std::pow(std::fabs(x[0]-h/2),alpha+beta)      // NOLINT
      //     * gsl_sf_hyperg_2F1(alpha+1, alpha+beta, alpha+beta+1, -x[j]/std::fabs(x[0]-h/2.0))    // NOLINT
      //     +                                                                                      // NOLINT
      //     y(x[jb])*s  td::pow(std::fabs(x[jb]),beta)/pow(std::fabs(x[n-1]+h/2),alpha+beta)       // NOLINT
      //     * gsl_sf_hyperg_2F1(alpha+1, alpha+beta, alpha+beta+1,  x[j]/std::fabs(x[n-1]+h/2.0))  // NOLINT
      //     );
      // }
  //    }

  void compute_far_field(std::vector<double> const & x ,
                         std::size_t                 ja,
                         std::size_t                 jb,
                         std::vector<double>       & ff) {  // NOLINT
    using boost::math::hypergeometric_pFq;

    std::size_t const n = x.size();
    double const C1a = b118::frlap::normal_const<1>(ealpha);
    double const ebeta = 1.0 - ealpha;
    double const h = x.at(2) - x.at(1);

    assert(ff.size() == n);

    for (std::size_t j = ja; j <= jb; ++j) {
      ff[j-ja] = -C1a/(ealpha+ebeta) * (
        y(x[ja]) * std::pow(std::fabs(x[ja]), ebeta)
                 / std::pow(std::fabs(x[0] - h/2), ealpha + ebeta)
                 * hypergeometric_pFq({ealpha + 1, ealpha + ebeta},
                                      {ealpha + ebeta + 1},
                                      - x[j] / std::fabs(x[0] - h/2))
        +
        y(x[jb]) * std::pow(std::fabs(x[jb]), ebeta)
                 / std::pow(std::fabs(x[n-1] + h/2), ealpha + ebeta)
                 * hypergeometric_pFq({ealpha + 1, ealpha + ebeta},
                                      {ealpha + ebeta + 1},
                                      + x[j] / std::fabs(x[n-1] + h/2)) );
      }
  }
}  // end namespace Benchmark701
