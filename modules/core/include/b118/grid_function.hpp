//   Copyright (C) 2024  Guilherme F. Fornel        <gffrnl@gmail.com>
//                       Fabio Souto de Azevedo     <fazedo@gmail.com>

#pragma once

#include <vector>
#include <stdexcept>
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

    
    // grid_function(grid_function const & other)
    //     : grid_function(other.thegrid)
    // {
    //     values = other.values;
    // }


    /*
    grid_function& operator=(grid_function const & other) {
        grid_function * tmp_ptr = this;
        grid_function * new_ptr = new grid_function(other.thegrid);
        new_ptr->values = other.values;

        this -> *new_ptr;
        new_ptr = tmp_ptr;
        delete new_ptr;
        
        return *this;
    }
    */

    // grid_function& operator=(grid_function && other) {
    //     grid_function f;
    //     f.thegrid = other.thegrid;
    //     f.values = other.values;
    //     return f;
    // }
    
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


namespace b118 { namespace distance {
template<typename Real>
Real max(grid_function<Real> const & f, grid_function<Real> const & g) {
    grid<Real> const x = f.get_grid();
    if (g.get_grid() != x)
        throw std::invalid_argument("g.get_grid != x");

    typename grid<Real>::size_type n = x.numnodes();

    double ret = 0.0;
    for (typename grid<Real>::size_type j = 0; j < n; ++j) {
        double const abs_dif = std::fabs(f[j] - g[j]);
        if (abs_dif > ret) ret = abs_dif;
    }
    return ret;
}
}}  // end namespace b118::distance
