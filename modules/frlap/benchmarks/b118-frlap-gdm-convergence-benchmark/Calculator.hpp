// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once

#include <iostream>

#include <memory>
#include <functional>
#include <chrono>
#include <stdexcept>
#include <b118/grid.hpp>
#include <b118/grid_function.hpp>
#include <b118/frlap/gdm/coefficients/method.hpp>
#include <b118/frlap/generalized_differences.hpp>

template<typename Real>
class Calculator {

 public:
    class Problem final : public std::enable_shared_from_this<Problem> {
        Problem(Real ealpha, b118::grid<double> grid,
                             b118::grid<double> viswin)
            : ealpha_(ealpha), grid_(grid), viswin_(viswin)
        {
            if (!viswin.is_subgrid(grid))
                throw std::invalid_argument("viswin is now a subgrid of grid");
        }

        Problem(Real ealpha, b118::grid<double> grid)
            : Problem(ealpha, grid, grid)
        {}

     public:
        virtual ~Problem() {};

        //template<template<typename> class Method, class F>
        template<class F>
        static std::shared_ptr<Problem> create(Real ealpha,
                                               b118::grid<double> grid,
                                               F y,
                                               b118::grid<double> viswin) {
            std::shared_ptr<Problem> prob_ptr
                    = std::make_shared<Problem>(Problem(ealpha, grid, viswin));
                      
            prob_ptr->y_ = b118::grid_function<Real>(prob_ptr->grid_, y);
            //prob_ptr->coeff_gtor_ = Method<Real>(ealpha, grid.spacing());
            
            return prob_ptr;
        }

        //template<template<typename> class Method, class F>
        template<class F>
        static std::shared_ptr<Problem> create(Real ealpha,
                                               b118::grid<double> grid,
                                               F y) {
            // return create<Method, F>(ealpha, grid, y, grid);
            return create(ealpha, grid, y, grid);
        }

        std::shared_ptr<Problem> getptr() {
            return this->shared_from_this();
        }

        Real ealpha() const { return ealpha_; }
        b118::grid<Real> grid() const { return grid_; }
        b118::grid<Real> viswin() const { return viswin_; }

        b118::grid_function<Real> const & y() const { return y_; }

     private:
        Real const ealpha_;
        b118::grid<Real> const grid_;
        b118::grid<Real> const viswin_;
        
        b118::grid_function<Real> y_;
    };


    
    // template<
    //     template<typename> class CoeffGenerator,
    //     template<class...> class FarFieldEstimatorKind,
    //     class ...Args
    //     >
    // void calculate(std::shared_ptr<Problem> const & prob_ptr,
    //                b118::grid_function<Real> * frlap_y_trunc,
    //                b118::grid_function<Real> * far_field_y,
    //                FarFieldEstimatorKind<Args...> ffkind) {
    //     using std::cout;
    //     using std::endl;

    //     // Computes the execution time
    //     auto time_0 = std::chrono::high_resolution_clock::now();
    //     //
    //     //  HERE GOES THE INITIALIZATION
    //     //
    //     //        b118::grid_function<Real> frlap_y_truncated(prob_ptr->viswin());
        
    //     b118::frlap::generalized_differences<
    //         Real,
    //         CoeffGenerator
    //         > method(prob_ptr->ealpha(), prob_ptr->grid());
        
    //     auto time_1 = std::chrono::high_resolution_clock::now();
    //     auto delta_1 =
    //             std::chrono::duration_cast<std::chrono::milliseconds>(time_1 - time_0);
    //     cout << "Duration of the initialization: "
    //          << delta_1.count() << "ms" << endl;
    //     //
    //     auto time_2 = std::chrono::high_resolution_clock::now();
    //     //
    //     //  HERE GOES THE COMPUTATION OF THE TRUNCATED FRLAP
    //     //
    //     method.compute_truncated(prob_ptr->y(), frlap_y_trunc);
    //     auto time_3 = std::chrono::high_resolution_clock::now();
    //     auto delta_2 =
    //             std::chrono::duration_cast<std::chrono::milliseconds>(time_3 - time_2);
    //     cout << "Duration of 'compute_truncated': " << delta_2.count() << "ms" << endl;
    //     //
    //     auto time_4 = std::chrono::high_resolution_clock::now();
    //     //
    //     //  HERE COMPUTE THE FAR-FIELD (IF NEEDED)
    //     //
    //     if (far_field_y != nullptr) {
    //         method.compute_far_field(ffkind, far_field_y);
    //     } else {
    //         method.compute_far_field(ffkind, frlap_y_trunc, true);
    //     }
    //     auto time_5 = std::chrono::high_resolution_clock::now();
    //     auto delta_3 =
    //             std::chrono::duration_cast<std::chrono::milliseconds>(time_5 - time_4);
    //     cout << "Duration of 'compute_far_field': " << delta_3.count() << "ms" << endl;
    //     //
    //     //  HERE WRITE THE RESULTS (IF NEEDED)
    //     //
    // }


    template<
        template<typename> class CoeffGenerator,
        template<class...> class FarFieldEstimatorKind,
        class ...Args
        >
    void calculate(std::shared_ptr<Problem> const & prob_ptr,
                   b118::grid_function<Real> * frlap_y_trunc,
                   b118::grid_function<Real> * far_field_y,
                   FarFieldEstimatorKind<Args...> ffkind) {
        using std::cout;
        using std::endl;

        // Computes the execution time
        auto time_0 = std::chrono::high_resolution_clock::now();
        //
        //  HERE GOES THE INITIALIZATION
        //
        //        b118::grid_function<Real> frlap_y_truncated(prob_ptr->viswin());
        
        b118::frlap::generalized_differences<
            Real,
            CoeffGenerator
            > method(prob_ptr->ealpha(), prob_ptr->grid());
        
        auto time_1 = std::chrono::high_resolution_clock::now();
        auto delta_1 =
                std::chrono::duration_cast<std::chrono::milliseconds>(time_1 - time_0);
        cout << "Duration of the initialization: "
             << delta_1.count() << "ms" << endl;
        //
        auto time_2 = std::chrono::high_resolution_clock::now();
        //
        //  HERE GOES THE COMPUTATION OF THE TRUNCATED FRLAP
        //
        method.compute_truncated(prob_ptr->y(), frlap_y_trunc);
        auto time_3 = std::chrono::high_resolution_clock::now();
        auto delta_2 =
                std::chrono::duration_cast<std::chrono::milliseconds>(time_3 - time_2);
        cout << "Duration of 'compute_truncated': " << delta_2.count() << "ms" << endl;
        //
        auto time_4 = std::chrono::high_resolution_clock::now();
        //
        //  HERE COMPUTE THE FAR-FIELD (IF NEEDED)
        //
        if (far_field_y != nullptr) {
            method.compute_far_field(ffkind, far_field_y);
        } else {
            method.compute_far_field(ffkind, frlap_y_trunc, true);
        }
        auto time_5 = std::chrono::high_resolution_clock::now();
        auto delta_3 =
                std::chrono::duration_cast<std::chrono::milliseconds>(time_5 - time_4);
        cout << "Duration of 'compute_far_field': " << delta_3.count() << "ms" << endl;
        //
        //  HERE WRITE THE RESULTS (IF NEEDED)
        //
    }

    // template<
    //     template<typename> class CoeffGenerator,
    //     template<class...> class FarFieldEstimatorKind,
    //     class ...Args
    //     >
    // void calculate(std::shared_ptr<Problem> const & prob_ptr,
    //                b118::grid_function<Real> * frlap_y,
    //                FarFieldEstimatorKind<Args...> ffkind) {
    //     calculate<CoeffGenerator>(prob_ptr, frlap_y, nullptr, ffkind);
    // }

    template<
        template<typename> class CoeffGenerator,
        template<class...> class FarFieldEstimatorKind,
        class ...Args
        >
    b118::grid_function<Real>
    frlap(std::shared_ptr<Problem> const & prob_ptr,
          FarFieldEstimatorKind<Args...> ffkind    ) {
        b118::grid_function<Real> frlap_y(prob_ptr->viswin());
        calculate<CoeffGenerator>(prob_ptr, &frlap_y, nullptr, ffkind);
        return frlap_y;
    }
};
