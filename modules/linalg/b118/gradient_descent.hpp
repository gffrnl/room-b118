#pragma once
#include <boost/math/tools/norms.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

namespace b118 { namespace linalg {

template<typename T>
std::size_t gradient_descent(
    ublas::symmetric_matrix<T, ublas::lower> const& A,
    ublas::vector<T>                         const& b,
    ublas::vector<T>                              & x,
    T tol = 10e-10) {
    
    using boost::math::tools::sup_norm;

    {
        //if (A.size1() != A.size2()) // this is guaranteed by symmetric_matrix
        //    throw ublas::bad_argument("A.size1() != A.size2()");
        if (b.size()  != A.size1())
            throw ublas::bad_argument("b.size() != A.size2()");
        if (x.size()  != A.size1())
            throw ublas::bad_argument("x.size() != A.size2()");
        if (!(tol > 0))
            throw ublas::bad_argument("tol <= 0");
    }

    ublas::vector<T> r(b.size());
    ublas::noalias(r) = b - ublas::prod(A, x);
    
    if (sup_norm(r) < tol) return 0;

    ublas::vector<T> Ar(b.size());
    std::size_t k = 1;
    for (;;) {
        T const rTr = ublas::inner_prod(r,r);
        ublas::noalias(Ar) = ublas::prod(A,r);
        T const gamma = rTr / ublas::inner_prod(r,Ar);
        x += gamma * r;
        if (rTr < tol) return k;
        r -= gamma * Ar;
        ++k;
    }

    return k;

} // end template gradient_descent()

}} // namespace b118::linalg