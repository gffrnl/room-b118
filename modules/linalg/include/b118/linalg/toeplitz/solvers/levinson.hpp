// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <cstddef>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <limits>

namespace b118 {
namespace linalg {
namespace toeplitz {
namespace solvers {

template<
    typename T,
    class InputBidirIt1, class InputBidirIt2,
    class OutputBidirIt
>
OutputBidirIt levinson(InputBidirIt1 row_beg, InputBidirIt1 row_end,
                       InputBidirIt2 rhs_beg,
                       OutputBidirIt sol_beg) {
    std::size_t n = 0;  // system's size

    {
        auto const len = std::distance(row_beg, row_end);

        if (len < static_cast<decltype(len)>(0))
            throw std::invalid_argument("levinson(): "
                                        "std::distance(row_beg, row_end) < 0");

        if (len == static_cast<decltype(len)>(0))
            return sol_beg;

        if (std::abs(*row_beg) < std::numeric_limits<T>::epsilon())
            throw std::invalid_argument("levinson(): zero diagonal");

        if (len == static_cast<decltype(len)>(1)) {
            *sol_beg = *rhs_beg / *row_beg;
            return std::next(sol_beg);
        }

        n = static_cast<std::size_t>(len);
    }

    std::vector<T> aux(n);  // auxiliary buffer

    auto row_it   = std::next(row_beg);
    auto row_rit  = --std::make_reverse_iterator(row_beg);
    auto rhs_it   = rhs_beg;
    auto sol_it   = sol_beg;
    auto aux_beg  = aux.begin();
    auto aux_it   = aux_beg;
    auto aux_rend = std::make_reverse_iterator(aux_beg);
    auto aux_rit  = aux_rend;

    T alpha, beta, mu;

    alpha = *aux_beg = - (*row_it / *row_beg);
    beta  = static_cast<T>(1);
    *sol_beg = (*rhs_beg / *row_beg);

    auto op_mu =
           [&mu](T const & a, T const & b) -> T { return a + mu * b;    };
    auto op_alpha =
        [&alpha](T const & a, T const & b) -> T { return alpha * a + b; };

    while (row_it != row_end) {
        ++row_it;
        --row_rit;
        ++rhs_it;
        ++sol_it;
        ++aux_it;
        --aux_rit;

        beta *= (static_cast<T>(1) - alpha * alpha);

        mu = (*rhs_it - std::transform_reduce(sol_beg, sol_it, row_rit,
                                              static_cast<T>(0)))
            / beta / *row_beg;

        std::transform(sol_beg, sol_it, aux_rit, sol_beg, op_mu);
        *sol_it = mu;

        if (row_it != row_end) {
            alpha = - (*row_it + std::transform_reduce(aux_beg, aux_it, row_rit,
                                                       static_cast<T>(0)))
                / beta / *row_beg;

            std::transform(aux_rit, aux_rend, aux_beg, std::prev(aux_rit),
                           op_alpha);
            *aux_beg = alpha;
            std::reverse(aux_beg, std::next(aux_it));
        }
    }

    return ++sol_it;
}

}  // end namespace solvers
}  // end namespace toeplitz
}  // end namespace linalg
}  // end namespace b118
