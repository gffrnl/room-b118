// Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once

#include <b118/linspace.hpp>
#include <b118/frlap/generalized_differences.hpp>
#include <b118/frlap/gdm/coefficients.hpp>
#include <b118/frlap/gdm/far_field_estimator.hpp>
#include <b118/grid.hpp>

#ifndef FRLAP_SIMPLE_CALCULATOR_FLOAT
typedef double real_t;
#else
typedef FRLAP_SIMPLE_CALCULATOR_FLOAT real_t
#endif

enum class CoefficientsKind {
    Spectral,
    SpectralQAWO,
    GorenfloMainardi,
    HuangOberman1,
    HuangOberman2,
    CPer3Point,
};

std::map<std::string, CoefficientsKind> const CoeffKindFromStr {
    {"spec"    , CoefficientsKind::Spectral        },
    {"specQAWO", CoefficientsKind::SpectralQAWO    },
    {"gormai"  , CoefficientsKind::GorenfloMainardi},
    {"huob1"   , CoefficientsKind::HuangOberman1   },
    {"huob2"   , CoefficientsKind::HuangOberman2   },
    {"cper3p"  , CoefficientsKind::CPer3Point      }
};

std::map<CoefficientsKind, std::string> const CoeffKindMessage {
    {CoefficientsKind::Spectral        ,
     "using `spectral` coefficients"         },
    {CoefficientsKind::SpectralQAWO    ,
     "using `spectral_qawo` coefficients"    },
    {CoefficientsKind::GorenfloMainardi,
     "using `gorenglo_mainardi` coefficients"},
    {CoefficientsKind::HuangOberman1   ,
     "using `huang_oberman_1` coefficients"  },
    {CoefficientsKind::HuangOberman2   ,
     "using `huang_oberman_2` coefficients"  },
    {CoefficientsKind::CPer3Point      ,
     "using `cper_3point` coefficients"      }
};


// enum class FarFieldKind {
//     Zero,
//     General,
//     Algebraic,
// };

// std::map<std::string, FarFieldKind> const FFKindFromStr {
//     {"zero"     , FarFieldKind::Zero     },
//     {"general"  , FarFieldKind::General  },
//     {"algebraic", FarFieldKind::Algebraic}
// };

// std::map<CoefficientsKind, std::string> const FFKindMessage {
//     {FarFieldKind::Zero     ,
//      "using zero far-field estimation"          },
//     {FarFieldKind::General  ,
//      "using general far-field estimator"        },
//     {FarFieldKind::Algebraic,
//      "using algebraic decay far-field estimator"}
// };

namespace b118::frlap::SimpleCalculator {
extern real_t ealpha; // frac. Laplacian exponent
extern real_t a0;     // left  inner endpoint
extern real_t b0;     // right inner endpoint
extern real_t y(real_t const & x);
extern real_t frLap_y(real_t const & x);
}  // end namespace b118::frlap::SimpleCalculator



