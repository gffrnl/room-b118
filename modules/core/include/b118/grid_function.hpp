//   Copyright (C) 2024  Guilherme F. Fornel        <gffrnl@gmail.com>
//                       Fabio Souto de Azevedo     <fazedo@gmail.com>

#pragma once

#include <vector>
#include "./grid.hpp"

namespace b118 {

template<typename Real = double>
class grid_function {
    grid<Real> thegrid;
    std::vector<Real> values;

 public:
    grid_function() {}

    template<class F>
    grid_function(grid<Real> const & g, F f)
        : thegrid(g), values(g.numnodes()) {
        typename grid<Real>::size_type n = thegrid.numnodes();
        for (typename grid<Real>::size_type k{0}; k < n; ++k)
#ifdef NDEBUG
            values[k] = f(thegrid[k]);
#else
            values.at(k) = f(thegrid[k]);
#endif
    }

    explicit grid_function(grid<Real> const & g)
        : thegrid(g), values(g.numnodes(), static_cast<Real>(0))
    {}

    inline Real operator[] (typename grid<Real>::size_type k) const {
#ifdef NDEBUG
        return values[k];
#else
        return values.at(k);
#endif
    }

    inline Real& operator[] (typename grid<Real>::size_type k) {
#ifdef NDEBUG
        return values[k];
#else
        return values.at(k);
#endif
    }

    inline auto data() const {
        return values.data();
    }

    inline auto data() {
        return values.data();
    }

    grid<Real> get_grid() const {
        return thegrid;
    }
};

}  // end namespace b118
