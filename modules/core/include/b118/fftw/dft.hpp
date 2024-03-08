// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once

#include "./plan.hpp"
#include "./complex_ptr.hpp"

namespace b118 {
namespace fftw {
namespace dft  {

namespace r2c {
    template<typename Real>
    fftw::plan<Real> create_plan(std::size_t n,
                                 Real * in,
                                 fftw::complex_ptr<Real> & out,  // NOLINT
                                 unsigned int flags) {
        fftw::plan<Real> p;
        p.raw = fftwq_plan_dft_r2c_1d(n, in, out.raw, flags);
        return p;
    }

    template<>
    fftw::plan<float> create_plan<float>(std::size_t n,
                                         float * in,
                                         fftw::complex_ptr<float> & out,  // NOLINT
                                         unsigned int flags) {
        fftw::plan<float> p;
        p.raw = fftwf_plan_dft_r2c_1d(n, in, out.raw, flags);
        return p;
    }

    template<>
    fftw::plan<double> create_plan<double>(std::size_t n,
                                           double * in,
                                           fftw::complex_ptr<double> & out,  // NOLINT
                                           unsigned int flags) {
        fftw::plan<double> p;
        p.raw = fftw_plan_dft_r2c_1d(n, in, out.raw, flags);
        return p;
    }

    template<>
    fftw::plan<long double> create_plan<long double>(std::size_t n,
                                           long double * in,
                                           fftw::complex_ptr<long double> & out,  // NOLINT
                                           unsigned int flags) {
        fftw::plan<long double> p;
        p.raw = fftwl_plan_dft_r2c_1d(n, in, out.raw, flags);
        return p;
    }


    template<typename Real>
    void execute(fftw::plan<Real> & p,  // NOLINT
                 Real * in,
                 fftw::complex_ptr<Real> & out) {  // NOLINT
        fftwq_execute_dft_r2c(p.raw, in, out.raw);
    }

    template<>
    void execute<float>(fftw::plan<float> & p,  // NOLINT
                        float * in,
                        fftw::complex_ptr<float> & out) {  // NOLINT
        fftwf_execute_dft_r2c(p.raw, in, out.raw);
    }

    template<>
    void execute<double>(fftw::plan<double> & p,  // NOLINT
                         double * in,
                         fftw::complex_ptr<double> & out) {  // NOLINT
        fftw_execute_dft_r2c(p.raw, in, out.raw);
    }

    template<>
    void execute<long double>(fftw::plan<long double> & p,  // NOLINT
                              long double * in,
                              fftw::complex_ptr<long double> & out) {  // NOLINT
        fftwl_execute_dft_r2c(p.raw, in, out.raw);
    }

}  // end namespace r2c


namespace c2r {
    template<typename Real>
    fftw::plan<Real> create_plan(std::size_t n,
                                 fftw::complex_ptr<Real> & in,   // NOLINT
                                 Real * out,
                                 unsigned int flags) {
        fftw::plan<Real> p;
        p.raw = fftwq_plan_dft_c2r_1d(n, in.raw, out, flags);
        return p;
    }

    template<>
    fftw::plan<float> create_plan<float>(std::size_t n,
                                         fftw::complex_ptr<float> & in,  // NOLINT
                                         float * out,
                                         unsigned int flags) {
        fftw::plan<float> p;
        p.raw = fftwf_plan_dft_c2r_1d(n, in.raw, out, flags);
        return p;
    }

    template<>
    fftw::plan<double> create_plan<double>(std::size_t n,
                                           fftw::complex_ptr<double> & in,  // NOLINT
                                           double * out,
                                           unsigned int flags) {
        fftw::plan<double> p;
        p.raw = fftw_plan_dft_c2r_1d(n, in.raw, out, flags);
        return p;
    }

    template<>
    fftw::plan<long double> create_plan<long double>(std::size_t n,
                                           fftw::complex_ptr<long double> & in,  // NOLINT
                                           long double * out,
                                           unsigned int flags) {
        fftw::plan<long double> p;
        p.raw = fftwl_plan_dft_c2r_1d(n, in.raw, out, flags);
        return p;
    }

    template<typename Real>
    void execute(fftw::plan<Real> & p,  // NOLINT
                 fftw::complex_ptr<Real> & in,  // NOLINT
                 Real * out) {
        fftwq_execute_dft_c2r(p.raw, in.raw, out);
    }

    template<>
    void execute<float>(fftw::plan<float> & p,  // NOLINT
                        fftw::complex_ptr<float> & in,  // NOLINT
                        float * out) {
        fftwf_execute_dft_c2r(p.raw, in.raw, out);
    }

    template<>
    void execute<double>(fftw::plan<double> & p,  // NOLINT
                        fftw::complex_ptr<double> & in, // NOLINT
                        double * out) {
        fftw_execute_dft_c2r(p.raw, in.raw, out);
    }

    template<>
    void execute<long double>(fftw::plan<long double> & p,  // NOLINT
                        fftw::complex_ptr<long double> & in, // NOLINT
                        long double * out) {
        fftwl_execute_dft_c2r(p.raw, in.raw, out);
    }
}  // end namespace c2r

}  // end namespace dft
}  // end namespace fftw
}  // end namespace b118
