// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>

#include <iostream>
#include "../include/b118/linalg/toeplitz/fsymmtp_manager.hpp"

int main() {
    using std::cout;
    using std::clog;
    using std::endl;
    using real_t = __float128;

    fsymmtp_manager<real_t>::initialize(11);
    fsymmtp_manager<real_t>::initialize(12);
    cout << "capacity = " << fsymmtp_manager<real_t>::capacity() << endl;
    cout << "size     = " << fsymmtp_manager<real_t>::size()     << endl;
    fsymmtp_manager<real_t>::initialize(13);
    cout << "capacity = " << fsymmtp_manager<real_t>::capacity() << endl;
    cout << "size     = " << fsymmtp_manager<real_t>::size()     << endl;

    fsymmtp_manager<real_t>::resize(16);
    cout << "capacity = " << fsymmtp_manager<real_t>::capacity() << endl;
    cout << "size     = " << fsymmtp_manager<real_t>::size()     << endl;

    fsymmtp_manager<real_t>::resize(32);
    cout << "capacity = " << fsymmtp_manager<real_t>::capacity() << endl;
    cout << "size     = " << fsymmtp_manager<real_t>::size()     << endl;

    return 0;
}