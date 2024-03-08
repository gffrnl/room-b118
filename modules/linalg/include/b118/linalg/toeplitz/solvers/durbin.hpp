// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <iterator>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace b118 {
namespace linalg {
namespace toeplitz {
namespace solvers {

template<
    typename T,
    class InputBidirIt,
    class OutputBidirIt
>
OutputBidirIt durbin(InputBidirIt  rhs_beg, InputBidirIt rhs_end,
                     OutputBidirIt sol_beg) {
    {
        auto const len = std::distance(rhs_beg, rhs_end);

        if (len < static_cast<decltype(len)>(0))
            throw std::invalid_argument("durbin(): "
                                        "std::distance(rhs_beg, rhs_end) < 0");

        if (len == static_cast<decltype(len)>(0))
            return sol_beg;

        if (len == static_cast<decltype(len)>(1)) {
            *sol_beg = *rhs_beg;
            return std::next(sol_beg);
        }
    }


    auto rhs_it = rhs_beg;
    auto sol_it = std::next(sol_beg);
    auto sol_rend = std::make_reverse_iterator(sol_beg);
    auto sol_rit  = sol_rend;

    T alpha, beta;

    alpha = *sol_beg = - (*rhs_it);
    beta = static_cast<T>(1);

    auto op =
        [&alpha](T const & a, T const & b) -> T { return alpha * a + b; };

    while (++rhs_it != rhs_end) {
        ++sol_it;
        --sol_rit;

        beta  *= (static_cast<T>(1) - alpha * alpha);

        alpha = - (*rhs_it + std::transform_reduce(rhs_beg, rhs_it, sol_rit,
                                                   static_cast<T>(0)))
            / beta;

        std::transform(sol_rit, sol_rend, sol_beg, std::prev(sol_rit),
                       op);
        *sol_beg = alpha;
        std::reverse(sol_beg, sol_it);
    }

    return ++sol_it;
}

}  // end namespace solvers
}  // end namespace toeplitz
}  // end namespace linalg
}  // end namespace b118
