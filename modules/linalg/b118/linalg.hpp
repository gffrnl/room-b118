#pragma once
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

namespace b118 { namespace linalg {

namespace matrix_kind {
    template<typename...> struct dense {};
    template<typename...> struct column {};
    template<typename...> struct symmetric_lower {};
    template<typename...> struct skew_symmetric {};
    template<typename...> struct sparse {};
    template<typename...> struct toeplitz {};
    template<typename...> struct symmetric_toeplitz {};
}


//template<typename T, class C>
//class matrix {};

template<typename T, template<typename...> class Kind = matrix_kind::dense>
class matrix;


template<typename T>
class matrix<T, matrix_kind::symmetric_lower>:
    public ublas::symmetric_matrix<T, ublas::lower> {

public: // exports
    using size_type       = typename ublas::symmetric_matrix<T, ublas::lower>::size_type;
    using value_type      = typename ublas::symmetric_matrix<T, ublas::lower>::value_type;
    using reference       = typename ublas::symmetric_matrix<T, ublas::lower>::reference;
    using const_reference = typename ublas::symmetric_matrix<T, ublas::lower>::const_reference;

public: // constructors
    matrix(size_type sz) : ublas::symmetric_matrix<T, ublas::lower>(sz) {}

public: // accessors
    inline const_reference operator()(size_type i, size_type j) const {
        return ublas::symmetric_matrix<T, ublas::lower>::operator()(i-1, j-1);
    }
    inline reference operator()(size_type i, size_type j) {
        return ublas::symmetric_matrix<T, ublas::lower>::operator()(i-1, j-1);
    }
};

template<typename T>
class matrix<T, matrix_kind::column>:
    public ublas::vector<T> {

public: // exports
    using size_type       = typename ublas::vector<T>::size_type;
    using value_type      = typename ublas::vector<T>::value_type;
    using reference       = typename ublas::vector<T>::reference;
    using const_reference = typename ublas::vector<T>::const_reference;

public: // constructors
    matrix(size_type sz) : ublas::vector<T>(sz) {}

public: // accessors
    inline const_reference operator()(size_type i) const {
        return ublas::vector<T>::operator()(i-1);
    }
    inline reference operator()(size_type i) {
        return ublas::vector<T>::operator()(i-1);
    }
    inline const_reference operator[](size_type i) const {
        return ublas::vector<T>::operator()(i);
    }
    inline reference operator[](size_type i) {
        return ublas::vector<T>::operator()(i);
    }
};

template<typename T>
using vector = matrix<T, matrix_kind::column>;


template<typename T>
class matrix<T, matrix_kind::symmetric_toeplitz> {
    std::vector<T> data;

public: // exports
    using size_type       = typename std::size_t;
    using value_type      = T;
    using reference       = T&;
    using const_reference = T const&;

public: // constructors
    matrix(size_type sz) : data(sz) {}

public: // destructor
    virtual ~matrix() {}

public: // accessors
    inline const_reference operator()(size_type i, size_type j) const {
#ifdef NDEBUG
        return data[(i > j)? i-j : j-i];
#else
        return data.at((i > j)? i-j : j-i);
#endif
    }
    inline reference operator()(size_type i, size_type j) {
#ifdef NDEBUG
        return data[(i > j)? i-j : j-i];
#else
        return data.at((i > j)? i-j : j-i);
#endif
    }

public:
    inline size_type size1() const { return data.size(); }
    inline size_type size2() const { return size1(); }
};

}}

template<class E, class T, class U>
std::basic_ostream<E, T>& operator<<(std::basic_ostream<E, T> &os,
                                     b118::linalg::matrix<U, b118::linalg::matrix_kind::symmetric_toeplitz> const& m) {
    using size_type = typename b118::linalg::matrix<U, b118::linalg::matrix_kind::symmetric_toeplitz>::size_type;
    size_type size1 = m.size1();
    size_type size2 = m.size2();
    std::basic_ostringstream<E, T, std::allocator<E> > s;
    s.flags(os.flags ());
    s.imbue(os.getloc ());
    s.precision(os.precision ());
    s << '[' << size1 << ',' << size2 << "](";
    if (size1 > 0) {
        s << '(' ;
        if (size2 > 0)
            s << m(0, 0);
        for (size_type j = 1; j < size2; ++j)
            s << ',' << m(0, j);
        s << ')';
    }
    for (size_type i = 1; i < size1; ++i) {
        s << ",(" ;
        if (size2 > 0)
            s << m(i, 0);
        for (size_type j = 1; j < size2; ++j)
            s << ',' << m(i, j);
        s << ')';
    }
    s << ')';
    return os << s.str().c_str();
}