// Copyright 2023  Guilherme F. Fornel
#include <cmath>
#include <iostream>
#include <b118/frlap.hpp>
#include <b118/frlap/gdm/trunc_uniform.hpp>
#include <b118/frlap/gdm/calculators/simple.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

using Interval = typename b118::frlap::gdm::calculators::Interval;
using Problem  = typename b118::frlap::gdm::calculators::Problem;

struct ProblemHuangOberman2014 : Problem {
    double y(double x) const override {
        double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
        return sgn * std::pow(1.0 + x*x, -0.5+ealpha/2.0);
    }

    double frLap_y(double x) const override {
        double const sgn = (ealpha > 1.0) ? -1.0 : 1.0;
        return sgn * std::exp2(ealpha) *
            std::tgamma(0.5+ealpha/2.0) / std::tgamma(0.5-ealpha/2.0) *
            pow(1.0 + x*x, -0.5-ealpha/2.0);
    }

    void compute_far_field(std::vector<double> const& x ,
                           std::size_t                ja,
                           std::size_t                jb,
                           double*                    ff) const override {
        using boost::math::hypergeometric_pFq;

        double const C1a = b118::frlap::normal_const<1>(ealpha);
        double const ebeta = 1.0 - ealpha;
        double const h = x[1] - x[0];

        for (std::size_t j = ja; j <= jb; ++j) {
            ff[j-ja] = -C1a/(ealpha+ebeta) * (
              y(x[ja]) * std::pow(std::fabs(x[ja]), ebeta)
                       / std::pow(std::fabs(x[0] - h/2), ealpha + ebeta)
                * hypergeometric_pFq({ealpha + 1, ealpha + ebeta},
                                     {ealpha + ebeta + 1},
                                     - x[j] / std::fabs(x[0] - h/2))
              +
              y(x[jb]) * std::pow(std::fabs(x[jb]), ebeta)
                       / std::pow(std::fabs(x[numnod-1] + h/2), ealpha + ebeta)
                * hypergeometric_pFq({ealpha + 1, ealpha + ebeta},
                                     {ealpha + ebeta + 1},
                                       x[j] / std::fabs(x[numnod-1] + h/2)) );
        }
    }
};

std::vector<std::unique_ptr<Problem>> make_prob_ptrs() {
    using b118::frlap::gdm_t;
    using b118::frlap::gdm::trunc_uniform_huob2;

    std::vector<std::unique_ptr<Problem>> prob_ptrs_vec;
    double const ealpha = 0.4;

    for (unsigned i = 0; i < 6; ++i) {
        for (unsigned j = 0; j < 4; ++j) {
            unsigned L0 = 4 << i;
            unsigned L  = 3 * L0;
            std::size_t numnod = (3 * 10 << (i+j)) + 1;
            std::string label;

            label = std::string("L").append(std::to_string(L0))
                                    .append("_")
                                    .append(std::to_string(numnod));

            std::unique_ptr<Problem> prob_ptr(new ProblemHuangOberman2014);

            prob_ptr->ealpha = ealpha;
            prob_ptr->domain = Interval{-static_cast<double>(L) ,
                            +static_cast<double>(L) };
            prob_ptr->viswin = Interval{-static_cast<double>(L0),
                            +static_cast<double>(L0)};
            prob_ptr->numnod = numnod;
            prob_ptr->method_ptr =
                std::unique_ptr<gdm_t>(new trunc_uniform_huob2);
            prob_ptr->label = label;

            prob_ptrs_vec.push_back(move(prob_ptr));
        }
    }

    prob_ptrs_vec.shrink_to_fit();
    return prob_ptrs_vec;
}


namespace b118  {
namespace frlap {
namespace gdm   {
namespace calculators  {

    std::vector<std::unique_ptr<Problem>> prob_ptrs = make_prob_ptrs();

}  // end namespace calculators
}  // end namespace gdm
}  // end namespace frlap
}  // end namespace b118
