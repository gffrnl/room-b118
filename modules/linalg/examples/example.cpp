/*
#include <b118/linalg.hpp>
#include <b118/linalg/conjugate_directions.hpp>
#include <b118/linalg/conjugate_gradient.hpp>
#include <b118/linalg/gradient_descent.hpp>
#include <b118/linalg/gram_schmidt.hpp>

int main() {
    using namespace std;
    using b118::linalg::matrix_kind::column;
    using b118::linalg::matrix_kind::dense;
    using b118::linalg::matrix_kind::sparse;
    using b118::linalg::matrix_kind::symmetric_lower;
    using b118::linalg::matrix_kind::symmetric_toeplitz;
    using b118::linalg::matrix_kind::toeplitz;

    // ublas::symmetric_matrix<double, ublas::lower> A(2);
    // ublas::vector<double> x(2), b(2);

    b118::linalg::matrix<double, symmetric_lower> A(2);
    // b118::linalg::matrix<double, column> x(2), b(2);
    b118::linalg::vector<double> x1(2), x2(2), b(2);

    A(1, 1) = 3.0;
    A(2, 1) = 2.0;
    A(2, 2) = 6.0;
    b(1) = 2.0;
    b(2) = -8.0;

    std::size_t niter;

    niter = b118::linalg::gradient_descent(A, b, x1);
    cout << "[\'gradient_descent\'] x = " << x1 << endl;
    cout << "\t [no. of iterations = " << niter << ']' << endl;

    niter = b118::linalg::conjugate_gradient(A, b, x2);
    cout << "[\'conjugate_gradient\'] x = " << x2 << endl;
    cout << "\t [no. of iterations = " << niter << ']' << endl;

    b118::linalg::matrix<double, symmetric_toeplitz> T(3);
    T(1, 1) = 10;
    T(2, 1) = 3;
    T(3, 1) = 1;
    T(2, 3) = 3.5;

    cout << "T = " << T << endl;

    return 0;
}
//*/

#include <iostream>
#include <complex>
#include <b118/linspace.hpp>
#include <b118/linalg.hpp>

int main() {
    using namespace std;
    vector<double> x;
    vector<complex<double>> xc;
    double step;
    complex<double> stepc;

    x = linspace<double, std::vector>(10.7, 44.1, 5);
    cout << "x = ";
    for (auto& elm : x) cout << elm << ", ";
    cout << endl;
    for (std::size_t k = 1; k < x.size(); ++k)
        cout << x[k] - x[k-1] << endl;
    cout << "------------" << endl;

    x = linspace<double, std::vector>(10.7, 44.1, 5, exclude::both);
    cout << "x = ";
    for (auto& elm : x) cout << elm << ", ";
    cout << endl;
    for (std::size_t k = 1; k < x.size(); ++k)
        cout << x[k] - x[k-1] << endl;
    cout << "------------" << endl;

    xc = linspace<std::complex<double>, std::vector>(0.0 + 0.0i, complex<double>(-1.0,1.0), 5, exclude::both);
    cout << "xc = ";
    for (auto& elm : xc) cout << elm << ", ";
    cout << endl;
    for (std::size_t k = 1; k < xc.size(); ++k)
        cout << xc[k] - xc[k-1] << endl;
    cout << "------------" << endl;


    step = linspace(x.begin(), x.end(), 1.1, 2.5, exclude::left);
    cout << "x = ";
    for (auto& elm : x) cout << elm << ", ";
    cout << endl;
    cout << "step = " << step << endl;
    for (std::size_t k = 1; k < x.size(); ++k)
        cout << x[k] - x[k-1] << endl;
    cout << "------------" << endl;


    step = linspace(x.rbegin(), x.rend(), 1.1, 2.5, exclude::left);
    cout << "x = ";
    for (auto& elm : x) cout << elm << ", ";
    cout << endl;
    cout << "step = " << step << endl;
    for (std::size_t k = 1; k < x.size(); ++k)
        cout << x[k] - x[k-1] << endl;
    cout << "------------" << endl;


    stepc = linspace(xc.begin(), xc.end(), 0.0 + 0.0i, complex<double>(-1.0,1.0), exclude::both);
    cout << "xc = ";
    for (auto& elm : xc) cout << elm << ", ";
    cout << endl;
    cout << "stepc = " << stepc << endl;
    for (std::size_t k = 1; k < xc.size(); ++k)
        cout << xc[k] - xc[k-1] << endl;
    cout << "------------" << endl;


    b118::linalg::vector<double> v(5);
    
    step = linspace<double>(v.begin(), v.end(), 1, 2.5, exclude::left | exclude::right);
    cout << "v = " << v << endl;
    cout << "step = " << step << endl;
    cout << "------------" << endl;

    step = linspace(v, 1.0, 2.5);
    cout << "v = " << v << endl;
    cout << "step = " << step << endl;
    cout << "------------" << endl;

    step = linspace(x, 1.0, 2.5);
    cout << "x = ";
    for (auto& elm : x) cout << elm << ", ";
    cout << endl;
    cout << "step = " << step << endl;
    for (std::size_t k = 1; k < x.size(); ++k)
        cout << x[k] - x[k-1] << endl;
    cout << "------------" << endl;


    auto u = linspace<double, b118::linalg::vector>(1.0, 2.5, 5, exclude::none, &step);
    cout << "u = " << u << endl;
    cout << "step = " << step << endl;
    cout << "------------" << endl;

    return 0;
}