// Copyright 2024 Guilherme F. Fornel <gffrnl@gmail.com>

#pragma once

#include <cstddef>
#include <cassert>
#ifdef FSYMTP_MANAGER_LOG
#include <iostream>
#endif

#include <b118/utility.hpp>
#include <b118/fftw/managed.hpp>
#include <b118/fftw/plan.hpp>

template<typename Real>
class fsymmtp_manager final : private b118::fftw::managed<Real> {
 public:
    static fsymmtp_manager& initialize(std::size_t n = 0) {
        static fsymmtp_manager obj(n);
        return obj;
    }

    static std::size_t capacity() {
        return initialize().m_capacity;
    }

    static std::size_t size() {
        return initialize().m_size;
    }

    static std::size_t padding() {
        return initialize().m_padding;
    }

    static std::size_t aug_size() {
        return initialize().m_aug_size;
    }

    static fsymmtp_manager & resize(std::size_t size) {
        return initialize().resize_impl(size);
    }

    static Real * const aug_row() {
        return initialize().m_aug_row;
    }

    static Real * const aug_vec() {
        return initialize().m_aug_vec;
    }

    static Real * const aux() {
        return initialize().m_aux;
    }

    static b118::fftw::plan<Real> const plan() {
        return initialize().m_plan;
    }

 private:
    fsymmtp_manager(fsymmtp_manager &) = delete;
    fsymmtp_manager& operator=(fsymmtp_manager&) = delete;

 private:
    std::size_t   m_capacity;
    std::size_t   m_size;
    std::size_t   m_padding;
    std::size_t   m_aug_size;

    Real * m_aug_row;
    Real * m_aug_vec;
    Real * m_aux;

    b118::fftw::plan<Real> m_plan;

    fsymmtp_manager(std::size_t n = 2) : m_capacity((n > 2)? n : 2),
                                         m_size(0),
                                         m_padding(0),
                                         m_aug_size(0),
                                         m_aug_row(nullptr),
                                         m_aug_vec(nullptr),
                                         m_aux(nullptr),
                                         m_plan() {
        if (n >= 2) {
            set_aug_size(m_capacity);
            alloc_fftw(m_aug_size);
        }
    }

    ~fsymmtp_manager() {
#ifdef FSYMTP_MANAGER_LOG
        std::clog << "*** fstp_manager destructor called ***"
                  << std::endl;
#endif
        destroy_plan_and_free();
    }

    void set_aug_size(std::size_t size) {
        m_padding  = (size % 2 == 0) ? (size - 4) / 2 : (size - 3) / 2;
        // (gffrnl) FFTW3 documentation says:
        //
        // https://www.fftw.org/fftw3_doc/
        //
        //
        // 4.3.5 Real-to-Real Transforms
        //
        // ...
        //
        // fftw_plan fftw_plan_r2r_1d(int n, double *in, double *out,
        //                            fftw_r2r_kind kind, unsigned flags);
        //
        // * n gives the (physical) size of the transform dimensions.
        //   They can be any positive integer.
        //
        //   - FFTW is generally best at handling sizes of the form
        //     2^a 3^b 5^c 7^d 11^e 13^f, where e+f is either 0 or 1, and the
        //     other exponents are arbitrary. Other sizes are computed by means
        //     of a slow, general-purpose algorithm (which nevertheless retains
        //     O(n log n) performance even for prime sizes). (It is possible to
        //     customize FFTW for different array sizes; see Installation and
        //     Customization.)
        //     Transforms whose sizes are powers of 2 are especially fast.
        //
        //   - For a REDFT00 or RODFT00 transform kind in a dimension of size n,
        //     it is n-1 or n+1, respectively, that should be factorizable in
        //     the above form.
        //
        // ...
        //
        // (gffrnl) So, since we use RODFT00, we choose the augmented size
        // (gffrnl) of the form 2^a - 1:
        m_aug_size = b118::next_exp2(size + 2 * (m_padding + 1)) - 1;
    }

    fsymmtp_manager & resize_impl(std::size_t size) {
        if (size != m_size) {
            std::size_t const old_aug_size = m_aug_size;
            set_aug_size(size);
            if (m_aug_size > old_aug_size) {
                destroy_plan_and_free();
                alloc_fftw(m_aug_size);
            }
            // resize the plan
            resize_plan();
            m_size = size;
            m_capacity = m_size;
        }
        return *this;
    }



    // To be defined to each type
    void alloc_fftw(std::size_t n);
    void resize_plan();
    void destroy_plan_and_free();
};


// FLOAT

template<>
void fsymmtp_manager<float>::destroy_plan_and_free() {
    if (m_plan.raw != nullptr) {
        fftwf_destroy_plan(m_plan.raw);
        m_plan.raw = nullptr;
    }
    fftwf_free(m_aux);
    m_aux = nullptr;
    fftwf_free(m_aug_vec);
    m_aug_vec = nullptr;
    fftwf_free(m_aug_row);
    m_aug_row = nullptr;
}

template<>
void fsymmtp_manager<float>::alloc_fftw(std::size_t n) {
    assert(m_aug_row == nullptr);
    assert(m_aug_vec == nullptr);
    assert(m_aux == nullptr);
    m_aug_row =
        static_cast<float *>(fftwf_malloc(n * sizeof(float)));
    m_aug_vec =
        static_cast<float *>(fftwf_malloc(n * sizeof(float)));
    m_aux =
        static_cast<float *>(fftwf_malloc(n * sizeof(float)));
    assert(m_aug_row != nullptr);
    assert(m_aug_vec != nullptr);
    assert(m_aux != nullptr);
}

template<>
void fsymmtp_manager<float>::resize_plan() {
    m_plan.raw = fftwf_plan_r2r_1d(m_aug_size, m_aug_row, m_aug_vec,
                                   FFTW_RODFT00, FFTW_ESTIMATE);
}


// DOUBLE

template<>
void fsymmtp_manager<double>::destroy_plan_and_free() {
    if (m_plan.raw != nullptr) {
        fftw_destroy_plan(m_plan.raw);
        m_plan.raw = nullptr;
    }
    fftw_free(m_aux);
    m_aux = nullptr;
    fftw_free(m_aug_vec);
    m_aug_vec = nullptr;
    fftw_free(m_aug_row);
    m_aug_row = nullptr;
}

template<>
void fsymmtp_manager<double>::alloc_fftw(std::size_t n) {
    assert(m_aug_row == nullptr);
    assert(m_aug_vec == nullptr);
    assert(m_aux == nullptr);
    m_aug_row =
        static_cast<double *>(fftw_malloc(n * sizeof(double)));
    m_aug_vec =
        static_cast<double *>(fftw_malloc(n * sizeof(double)));
    m_aux =
        static_cast<double *>(fftw_malloc(n * sizeof(double)));
    assert(m_aug_row != nullptr);
    assert(m_aug_vec != nullptr);
    assert(m_aux != nullptr);
}

template<>
void fsymmtp_manager<double>::resize_plan() {
    m_plan.raw = fftw_plan_r2r_1d(m_aug_size, m_aug_row, m_aug_vec,
                                  FFTW_RODFT00, FFTW_ESTIMATE);
}


// LONG DOUBLE

template<>
void fsymmtp_manager<long double>::destroy_plan_and_free() {
    if (m_plan.raw != nullptr) {
        fftwl_destroy_plan(m_plan.raw);
        m_plan.raw = nullptr;
    }
    fftwl_free(m_aux);
    m_aux = nullptr;
    fftwl_free(m_aug_vec);
    m_aug_vec = nullptr;
    fftwl_free(m_aug_row);
    m_aug_row = nullptr;
}

template<>
void fsymmtp_manager<long double>::alloc_fftw(std::size_t n) {
    assert(m_aug_row == nullptr);
    assert(m_aug_vec == nullptr);
    assert(m_aux == nullptr);
    m_aug_row =
        static_cast<long double *>(fftwl_malloc(n * sizeof(long double)));
    m_aug_vec =
        static_cast<long double *>(fftwl_malloc(n * sizeof(long double)));
    m_aux =
        static_cast<long double *>(fftwl_malloc(n * sizeof(long double)));
    assert(m_aug_row != nullptr);
    assert(m_aug_vec != nullptr);
    assert(m_aux != nullptr);
}

template<>
void fsymmtp_manager<long double>::resize_plan() {
    m_plan.raw = fftwl_plan_r2r_1d(m_aug_size, m_aug_row, m_aug_vec,
                                   FFTW_RODFT00, FFTW_ESTIMATE);
}


// DEFAULT - QUADRUPLE (__float128)

template<typename Real>
void fsymmtp_manager<Real>::destroy_plan_and_free() {
    if (m_plan.raw != nullptr) {
        fftwq_destroy_plan(m_plan.raw);
        m_plan.raw = nullptr;
    }
    fftwq_free(m_aux);
    m_aux = nullptr;
    fftwq_free(m_aug_vec);
    m_aug_vec = nullptr;
    fftwq_free(m_aug_row);
    m_aug_row = nullptr;
}

template<typename Real>
void fsymmtp_manager<Real>::alloc_fftw(std::size_t n) {
    assert(m_aug_row == nullptr);
    assert(m_aug_vec == nullptr);
    assert(m_aux == nullptr);
    m_aug_row =
        static_cast<Real *>(fftwq_malloc(n * sizeof(Real)));
    m_aug_vec =
        static_cast<Real *>(fftwq_malloc(n * sizeof(Real)));
    m_aux =
        static_cast<Real *>(fftwq_malloc(n * sizeof(Real)));
    assert(m_aug_row != nullptr);
    assert(m_aug_vec != nullptr);
    assert(m_aux != nullptr);
}

template<typename Real>
void fsymmtp_manager<Real>::resize_plan() {
    m_plan.raw = fftwq_plan_r2r_1d(m_aug_size, m_aug_row, m_aug_vec,
                                   FFTW_RODFT00, FFTW_ESTIMATE);
}
