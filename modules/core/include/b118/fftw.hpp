/*   libb118
 *
 *   modules/core/include/b118/fftw.hpp
 *
 *   C++ wrappers for fftw3
 *
 *   Copyright (C) 2024  Guilherme F. Fornel <gffrnl@gmail.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <fftw3.h>
#include <cstddef>
#ifdef B118_FFTW_MANAGER_LOG
    #include <iostream>
#endif

//  fftw3.h: line 431
//  __float128 (quad precision) is a gcc extension on i386, x86_64, and ia64
//  for gcc >= 4.6 (compiled in FFTW with --enable-quad-precision)
#if  (  __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)         )    && \
    !(  defined(__ICC)      || defined(__INTEL_COMPILER) || \
        defined(__CUDACC__) || defined(__PGI)                          )    && \
     (  defined(__i386__) || defined(__x86_64__) || defined(__ia64__)  )
    #define B118_CORE_FFTW3_ENABLE_QUAD
#endif


namespace b118 { namespace fftw {

    namespace detail {

        template<typename Real>
        struct complex_ctype;

        template<>
        struct complex_ctype<float>       final { using type = fftwf_complex; };

        template<>
        struct complex_ctype<double>      final { using type = fftw_complex;  };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
        template<typename Real>
        struct complex_ctype              final { using type = fftwl_complex; };
#else
        template<>
        struct complex_ctype<long double> final { using type = fftwl_complex; };
#endif

#if  defined(B118_CORE_FFTW3_ENABLE_QUAD)
        template<typename Real>
        struct complex_ctype              final { using type = fftwq_complex; };
#endif

    }  // end namespace detail


    // A C++-style syntatic sugar for fftw complex
    template<typename Real>
    using complex = typename detail::complex_ctype<Real>::type;


    // Wrappers for fftw3 malloc
    template<typename Real>
    struct alloc;

    template<>
    struct alloc<float> final {
        static float * real(std::size_t n) {
            return
                static_cast<float *>(
                    fftwf_malloc(static_cast<size_t>(n) * sizeof(float))
                );
        }
        static fftw::complex<float> * complex(std::size_t n) {
            return
                static_cast<fftwf_complex *>(
                    fftwf_malloc(static_cast<size_t>(n) * sizeof(fftwf_complex))
                );
        }
    };

    template<>
    struct alloc<double> final {
        static double * real(std::size_t n) {
            return
                static_cast<double *>(
                    fftw_malloc(static_cast<size_t>(n) * sizeof(double))
                );
        }
        static fftw::complex<double> * complex(std::size_t n) {
            return
                static_cast<fftw_complex *>(
                    fftw_malloc(static_cast<size_t>(n) * sizeof(fftw_complex))
                );
        }
    };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    struct alloc final
#else
    template<>
    struct alloc<long double> final
#endif
    {
        static long double * real(std::size_t n) {
            return
                static_cast<long double *>(
                    fftwl_malloc(static_cast<size_t>(n) * sizeof(long double))
                );
        }
        static fftw::complex<long double> * complex(std::size_t n) {
            return
                static_cast<fftwl_complex *>(
                    fftwl_malloc(static_cast<size_t>(n) * sizeof(fftwl_complex))
                );
        }
    };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    struct alloc final {
        static __float128 * real(std::size_t n) {
            return
                static_cast<__float128 *>(
                    fftwq_malloc(static_cast<size_t>(n) * sizeof(__float128))
                );
        }
        static fftw::complex<__float128> * complex(std::size_t n) {
            return
                static_cast<fftwq_complex *>(
                    fftwq_malloc(static_cast<size_t>(n) * sizeof(fftwq_complex))
                );
        }
    };
#endif


    // Wrappers for fftw3 free
    template<typename Real>
    struct free;

    template<>
    struct free<float> final {
        static void real(float * p) {
            fftwf_free(p);
        }
        static void complex(fftw::complex<float> * p) {
            fftwf_free(p);
        }
    };

    template<>
    struct free<double> final {
        static void real(double * p) {
            fftw_free(p);
        }
        static void complex(fftw::complex<double> * p) {
            fftw_free(p);
        }
    };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    struct free final
#else
    template<>
    struct free<long double> final
#endif
    {
        static void real(long double * p) {
            fftwl_free(p);
        }
        static void complex(fftw::complex<long double> * p) {
            fftwl_free(p);
        }
    };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    struct free final {
        static void real(__float128 * p) {
            fftwq_free(p);
        }
        static void complex(fftw::complex<__float128> * p) {
            fftwq_free(p);
        }
    };
#endif


    // fftw3 plans

    namespace detail {

        template<typename Real>
        struct plan_destroyer;


        template<>
        struct plan_destroyer<float> {
            static void destroy(fftwf_plan pl) {
                if (pl != NULL) {
                    fftwf_destroy_plan(pl);
                    pl = NULL;
                }
            }
        };

        template<>
        struct plan_destroyer<double> {
            static void destroy(fftw_plan pl) {
                if (pl != NULL) {
                    fftw_destroy_plan(pl);
                    pl = NULL;
                }
            }
        };

        template<>
        struct plan_destroyer<long double> {
            static void destroy(fftwl_plan pl) {
                if (pl != NULL) {
                    fftwl_destroy_plan(pl);
                    pl = NULL;
                }
            }
        };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
        template<>
        struct plan_destroyer<__float128> {
            static void destroy(fftwq_plan pl) {
                if (pl != NULL) {
                    fftwq_destroy_plan(pl);
                    pl = NULL;
                }
            }
        };
#endif


    }  // end namespace detail


    enum plan_kind { r2r, r2c, c2r };
    
    enum class r2rkind {
        redft00 = FFTW_REDFT00,
        redft01 = FFTW_REDFT01,
        redft10 = FFTW_REDFT10,
        redft11 = FFTW_REDFT11,
        rodft00 = FFTW_RODFT00,
        rodft01 = FFTW_RODFT01,
        rodft10 = FFTW_RODFT10,
        rodft11 = FFTW_RODFT11,
    };

    // Wrappers for fftw3 plans
    template<typename Real, unsigned Kind>
    class plan;


    // r2r plans
    // =========

    template<>
    class plan<float, r2r> final {
        fftwf_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             float * in, float * out,
             enum class r2rkind kind,
             unsigned int flags) {
            raw = fftwf_plan_r2r_1d(static_cast<size_t>(n),
                                    in, out,
                                    static_cast<fftw_r2r_kind>(kind),
                                    flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<float>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<float>::destroy(raw);
        }

        void execute(float * in, float * out) const {
            fftwf_execute_r2r(raw, in, out);
        }
    };

    template<>
    class plan<double, r2r> final {
        fftw_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             double * in, double * out,
             enum class r2rkind kind,
             unsigned int flags) {
            raw = fftw_plan_r2r_1d(static_cast<size_t>(n),
                                   in, out,
                                   static_cast<fftw_r2r_kind>(kind),
                                   flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<double>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }

        // Destructor
        ~plan() {
            detail::plan_destroyer<double>::destroy(raw);
        }

        void execute(double * in, double * out) const {
            fftw_execute_r2r(raw, in, out);
        }
    };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    class plan<Real, r2r> final
#else
    template<>
    class plan<long double, r2r> final
#endif
    {
        fftwl_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             long double * in, long double * out,
             enum class r2rkind kind,
             unsigned int flags) {
            raw = fftwl_plan_r2r_1d(static_cast<size_t>(n),
                                    in, out,
                                    static_cast<fftw_r2r_kind>(kind),
                                    flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<long double>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<long double>::destroy(raw);
        }

        void execute(long double * in, long double * out) const {
            fftwl_execute_r2r(raw, in, out);
        }
    };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<>
    class plan<__float128, r2r> final {
        fftwq_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             __float128 * in, __float128 * out,
             enum class r2rkind kind,
             unsigned int flags) {
            raw = fftwq_plan_r2r_1d(static_cast<size_t>(n),
                                    in, out,
                                    static_cast<fftw_r2r_kind>(kind),
                                    flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<__float128>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<__float128>::destroy(raw);
        }

        void execute(__float128 * in, __float128 * out) const {
            fftwq_execute_r2r(raw, in, out);
        }
    };
#endif


    // r2c plans
    // =========

    template<>
    class plan<float, r2c> final {
        fftwf_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             float * in, fftw::complex<float> * out,
             unsigned int flags) {
            raw = fftwf_plan_dft_r2c_1d(static_cast<size_t>(n),
                                        in, out,
                                        flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<float>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }

        // Destructor
        ~plan() {
            detail::plan_destroyer<float>::destroy(raw);
        }

        void execute(float * in, fftw::complex<float> * out) const {
            fftwf_execute_dft_r2c(raw, in, out);
        }
    };

    template<>
    class plan<double, r2c> final {
        fftw_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             double * in, fftw::complex<double> * out,
             unsigned int flags) {
            raw = fftw_plan_dft_r2c_1d(static_cast<size_t>(n),
                                       in, out,
                                       flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<double>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<double>::destroy(raw);
        }

        void execute(double * in, fftw::complex<double> * out) const {
            fftw_execute_dft_r2c(raw, in, out);
        }
    };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    class plan<Real, r2c> final
#else
    template<>
    class plan<long double, r2c> final
#endif
    {
        fftwl_plan raw;
    
     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             long double * in, fftw::complex<long double> * out,
             unsigned int flags) {
            raw = fftwl_plan_dft_r2c_1d(static_cast<size_t>(n),
                                        in, out,
                                        flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<long double>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<long double>::destroy(raw);
        }

        void execute(long double * in, fftw::complex<long double> * out) const {
            fftwl_execute_dft_r2c(raw, in, out);
        }
    };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    class plan<Real, r2c> final {
        fftwq_plan raw;

     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             __float128 * in, fftw::complex<__float128> * out,
             unsigned int flags) {
            raw = fftwq_plan_dft_r2c_1d(static_cast<size_t>(n),
                                        in, out,
                                        flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<__float128>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<__float128>::destroy(raw);
        }

        void execute(__float128 * in, fftw::complex<__float128> * out) const {
            fftwq_execute_dft_r2c(raw, in, out);
        }
    };
#endif

    // c2r plans
    // =========

    template<>
    class plan<float, c2r> final {
        fftwf_plan raw;

     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             fftw::complex<float> * in, float * out,
             unsigned int flags) {
            raw = fftwf_plan_dft_c2r_1d(static_cast<size_t>(n),
                                        in, out,
                                        flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<float>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<float>::destroy(raw);
        }

        void execute(fftw::complex<float> * in, float * out) const {
            fftwf_execute_dft_c2r(raw, in, out);
        }
    };

    template<>
    class plan<double, c2r> final {
        fftw_plan raw;

     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             fftw::complex<double> * in, double * out,
             unsigned int flags) {
            raw = fftw_plan_dft_c2r_1d(static_cast<size_t>(n),
                                       in, out,
                                       flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<double>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<double>::destroy(raw);
        }

        void execute(fftw::complex<double> * in, double * out) const {
            fftw_execute_dft_c2r(raw, in, out);
        }
    };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    class plan<Real, c2r> final
#else
    template<>
    class plan<long double, c2r> final
#endif
    {
        fftwl_plan raw;

     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             fftw::complex<long double> * in, long double * out,
             unsigned int flags) {
            raw = fftwl_plan_dft_c2r_1d(static_cast<size_t>(n),
                                        in, out,
                                        flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<long double>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<long double>::destroy(raw);
        }

        void execute(fftw::complex<long double> * in, long double * out) const {
            fftwl_execute_dft_c2r(raw, in, out);
        }
    };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    class plan<Real, c2r> final {
        fftwq_plan raw;

     public:
        // Deleted copy constructor and copy assignment:
        plan(plan const &) = delete;
        plan& operator=(plan const &) = delete;

        // Default constructor
        plan() { raw = NULL; }

        // Parametrized constructor
        plan(std::size_t n,
             fftw::complex<__float128> * in, __float128 * out,
             unsigned int flags) {
            raw = fftwq_plan_dft_c2r_1d(static_cast<size_t>(n),
                                        in, out,
                                        flags);
        }

        // Move constructor
        plan(plan && other) {
            raw = other.raw;
            other.raw = NULL;
        }

        // Move assignment
        plan& operator=(plan && other) {
            detail::plan_destroyer<__float128>::destroy(raw);
            raw = other.raw;
            other.raw = NULL;
            return *this;
        }
        
        // Destructor
        ~plan() {
            detail::plan_destroyer<__float128>::destroy(raw);
        }

        void execute(fftw::complex<__float128> * in, __float128 * out) const {
            fftwq_execute_dft_c2r(raw, in, out);
        }
    };
#endif


    // A fftw manager to cleanup system information

    template<typename Real>
    struct manager;

    template<>
    struct manager<float> final {
        static void cleanup() {
            fftwf_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
            std::clog << "*** fftwf_cleanup() called ***" << std::endl;
#endif
        }
    };

    template<>
    struct manager<double> final {
        static void cleanup() {
            fftw_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
            std::clog << "*** fftw_cleanup() called ***" << std::endl;
#endif
        }
    };

#if !defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    struct manager final
#else
    template<>
    struct manager<long double> final
#endif
    {
        static void cleanup() {
            fftwl_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
            std::clog << "*** fftwl_cleanup() called ***" << std::endl;
#endif
        }
    };

#if defined(B118_CORE_FFTW3_ENABLE_QUAD)
    template<typename Real>
    struct manager final {
        static void cleanup() {
            fftwq_cleanup();
#ifdef B118_FFTW_MANAGER_LOG
            std::clog << "*** fftwq_cleanup() called ***" << std::endl;
#endif
        }
    };
#endif


    // A fftw managed class whose fftw3-managed classes MUST inherit

    // template<typename Real>
    // class managed {
    //     unsigned num_users;

    // protected:
    //     managed() {
    //         static unsigned instances = 0;
    //         num_users = ++instances;
    //     }

    //     virtual ~managed() {
    //         if (--num_users == 0) b118::fftw::manager<Real>::cleanup();
    //         std::cout << "fftw_managed users = " << num_users << std::endl;
    //     }

    // public:
    //     unsigned users() const { return num_users; }
    // };

    template<typename Real>
    class managed {
        static unsigned instances;

    protected:
        managed() {
            ++instances;
        }

        virtual ~managed() {
            if (--instances == 0)
                manager<Real>::cleanup();
        }

    public:
        static unsigned num_users() { return instances; }
    };

    template<typename Real>
    unsigned managed<Real>::instances = 0;

}}  // end namespace b118::fftw
