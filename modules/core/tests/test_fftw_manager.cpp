// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>

#include <iostream>
#define FFTW_MANAGER_LOG
#include "../include/b118/fftw_manager.hpp"

template<typename Real>
struct Foo final : private fftw_managed<Real> {
    // using fftw_managed<Real>::users;

    ~Foo() {
        std::cout << "*** ~Foo() called ***" << std::endl;
    }

    inline unsigned fftw_managed_users() const {
        return fftw_managed<Real>::users();
    }
};

template<typename Real>
void bar() {
    using std::cout;
    using std::clog;
    using std::endl;

    Foo<Real> foo;
    cout << "foo.fftw_managed_users() " << foo.fftw_managed_users() << endl;
}

int main() {
    using std::cout;
    using std::clog;
    using std::endl;
    using real_t = __float128;

    clog << "\t***test_fftw_manager ***\n" << endl;

    Foo<real_t> foo1;
    Foo<real_t> foo2;
    {
        Foo<real_t> foo3;
        cout << "foo2.fftw_managed_users() "
             << foo2.fftw_managed_users()
             << endl;
        cout << "foo3.fftw_managed_users() "
             << foo3.fftw_managed_users()
             << endl;
    }
    bar<real_t>();
    cout << "foo2.fftw_managed_users() "
         << foo2.fftw_managed_users()
         << endl;

    return 0;
}
