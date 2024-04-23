// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once
#include <variant>

namespace ConvergenceBenchmark {
    enum class Method {
        spectral,
        spectral_qawo,
        gorenflo_mainardi,
        huang_oberman_1,
        huang_oberman_2,
        cper3point
    };

    namespace FarField {
        constexpr unsigned zero = 0;

        constexpr unsigned general = 1;

        template<typename Real>
        struct algebraic final {
            static constexpr unsigned id = 2;
            operator unsigned() const { return id; }
            Real const decay;
            explicit algebraic(Real decay) : decay(decay) {}
        };
    };

    template<typename Real>
    using method_t   = Method;

    template<typename Real>
    using farfield_t = std::variant<unsigned, FarField::algebraic<Real>>;

    template<typename Real>
    struct Problem {
        Real             const alpha;
        method_t<Real>   const method;
        farfield_t<Real> const farfield;

        Problem(Real alpha, method_t<Real> method, farfield_t<Real> farfield)
            : alpha(alpha), method(method), farfield(farfield)
        {}

        virtual ~Problem() {}
        virtual Real y(double const & x) = 0;
        virtual Real frLap_y(double const & x) = 0;
   
    public:
        std::string get_method_label() {
            std::string label;
            
            switch (method) {
                case (Method::spectral):
                    label = std::string("spectral");
                    break;

                case (Method::spectral_qawo):
                    label = std::string("spectral_qawo");
                    break;

                case (Method::gorenflo_mainardi):
                    label = std::string("gorenflo_mainardi");
                    break;

                case (Method::huang_oberman_1):
                    label = std::string("huang_oberman_1");
                    break;

                case (Method::huang_oberman_2):
                    label = std::string("huang_oberman_2");
                    break;

                case (Method::cper3point):
                    label = std::string("cper3point");
                    break;
                    
                default:
                    label = std::string("unknown");
            }

            return label;
        }

        std::string get_farfield_label() {
            std::string label;

            if (std::holds_alternative<FarField::algebraic<Real>>(farfield)) {
                label = std::string("algebraic");
            } else {
                switch (std::get<unsigned>(farfield)) {
                    case FarField::zero:
                        label = std::string("zero");
                        break;
                    case FarField::general:
                        label = std::string("general");
                        break;
                    default:
                        label = std::string("unknown");
                }
            }
            
            return label;
        }

        std::string get_label() {
            return (get_method_label() + "_" + get_farfield_label());
        }
        
        Real get_decay() const {
            if (std::holds_alternative<FarField::algebraic<Real>>(farfield))
                return std::get<FarField::algebraic<Real>>(farfield).decay;

            return -1;
        }
    };


    template<typename Real>
    using Problem_ptr = std::shared_ptr<Problem<Real>>;

    template<typename Real>
    Problem_ptr<Real> create_problem();

}  // end namespace ConvergenceBenchmark
