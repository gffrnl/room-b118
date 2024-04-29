//   Copyright (C) 2024  Guilherme F. Fornel        <gffrnl@gmail.com>
//                       Fabio Souto de Azevedo     <fazedo@gmail.com>

#pragma once

#include <cstddef>
#include <cassert>  // TODO(gffrnl): change to exceptions;
#include <utility>
#include <vector>
#include "./b118/almost_equal.hpp"


namespace b118 {

template<typename Real = double>
class grid {
 public:
    using size_type = std::size_t;

 private:
    Real const a;
    Real const b;
    size_type const n;
    Real const step;

 public:
    struct iterator {
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type        = Real;
        using difference_type   = std::ptrdiff_t;
        using pointer           = void;
        using reference         = Real &;

        grid const * const that;
        size_type n;

     public:
        iterator& operator++() {
            ++n;
            return *this;
        }

        iterator operator++(int) {
            iterator ret = *this;
            ++(*this);
            return ret;
        }

        iterator& operator--() {
            --n;
            return *this;
        }

        iterator operator--(int) {
            iterator ret = *this;
            --(*this);
            return ret;
        }

        Real operator*() const {
            return (*that)[n];
        }

        bool operator!=(iterator const& rhs) const {
            return n != rhs.n;
        }

        Real operator[] (size_type k) const {
            return (*that)[n + k];
        }
    };

 public:
    grid() : a(0), b(0), n(0), step(0) {}

    grid(std::pair<Real, Real> endpoints, size_type numnodes)
        : a(endpoints.first), b(endpoints.second), n(numnodes),
          step((b - a)/(n - 1))
    {}

    /*
    grid& operator=(grid const & other) {
        grid g(std::make_pair(other.a, other.b), other.n);
        std::swap(g, *this);
        return *this;
    }
    */
    
    inline Real lendpoint() const { return a; }
    inline Real rendpoint() const { return b; }
    
    inline size_type numnodes() const { return n; }

    inline Real spacing() const { return step; }

    inline Real operator[] (size_type k) const {
        Real const x = a / static_cast<Real>(n - 1);
        Real const y = b / static_cast<Real>(n - 1);
        return x * static_cast<Real>(n - 1 - k) +  y * static_cast<Real>(k);
    }

    bool operator== (grid const & other) const {
        if (other.a != a) return false;
        if (other.b != b) return false;
        if (other.n != n) return false;
        return true;
    }

    bool operator!= (grid const & other) const {
        return !(*this == other);
    }

    iterator begin() const {
         return iterator{this, static_cast<size_type>(0)};
    }

    iterator end() const {
        return iterator{this, n};
    }

    size_type closest(Real const & x) const {
        size_type pos =
            static_cast<size_type>((x - a) / step + static_cast<Real>(1)/2);
        if (pos < 0)
            return 0;
        else if (pos >= n)
            return n;
        return pos;
    }

    grid<Real> subgrid(size_type ka, size_type kb) {
        // preconditions:
        assert(ka >= 0);
        assert(kb <= n - 1);
        assert(ka <= kb);

        return grid<Real>({(*this)[ka], (*this)[kb]}, kb - ka + 1);
    }

    grid<Real> subgrid(Real a0, Real b0) {
        // preconditions:
        assert(!(a0 < a));
        assert(!(b0 > b));
        assert(!(a0 > b0));

        return subgrid(this->closest(a0), this->closest(b0));
    }

    inline grid<Real> subgrid(std::pair<Real, Real> endpoints) {
        return subgrid(endpoints.first, endpoints.second);
    }

    bool has_subgrid(grid<Real> const & other) const {
        if (!b118::almost_equal(other.step, step))
            return false;

        size_type const ka = (*this).closest(other[0]);
        if (!b118::almost_equal((*this)[ka], other[0]))
            return false;

        size_type const kb = (*this).closest(other[other.numnodes()-1]);
        if ((*this)[kb] != other[other.numnodes()-1])
            return false;

        return true;
    }

    bool is_subgrid(grid<Real> const & other) const {
        return other.has_subgrid(*this);
    }
};

}  // end namespace b118
