// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

//#include "./Calculator.hpp"

#include <cmath>
#include <memory>
#include <iostream>
#include "./ConvergenceBenchmark.hpp"

struct MyProblem final : public ConvergenceBenchmark::Problem<double> {
    MyProblem(double                                   alpha   ,
              ConvergenceBenchmark::method_t<double>   method  ,
              ConvergenceBenchmark::farfield_t<double> farfield)
        : Problem(alpha, method, farfield),
          sgn((alpha > 1.0) ? -1.0 : 1.0)
    {}

    double y(double const & x) override {   
        return sgn * pow(1.0 + x*x, -0.5 + alpha/2.0);
    }

    double frLap_y(double const & x) override {
        return sgn * pow(1.0 + x*x, -0.5 - alpha/2.0) * exp2(alpha)
                   * tgamma(0.5+alpha/2.0) / tgamma(0.5-alpha/2.0);
    }
    
 private:
    double const sgn;
};


namespace ConvergenceBenchmark {
    template<>
    Problem_ptr<double> create_problem<double>() {
        double const alpha = 0.4;
        double const beta  = 1.0 - alpha;

        return
                std::make_shared<MyProblem>(
                    alpha,
                    //ConvergenceBenchmark::Method::spectral,
                    //ConvergenceBenchmark::Method::spectral_qawo,
                    //ConvergenceBenchmark::Method::gorenflo_mainardi,
                    //ConvergenceBenchmark::Method::huang_oberman_1,
                    ConvergenceBenchmark::Method::huang_oberman_2,
                    //ConvergenceBenchmark::Method::cper3point,
                    //
                    //ConvergenceBenchmark::FarField::zero
                    ConvergenceBenchmark::FarField::general
                    //ConvergenceBenchmark::FarField::algebraic(beta)
                );
    }
}  // end namespace ConvergenceBenchmark
