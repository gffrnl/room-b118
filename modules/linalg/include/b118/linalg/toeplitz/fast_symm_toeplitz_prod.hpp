/*   libb118
 *
 *   modules/linalg/include/b118/linalg/toeplitz/fast_symm_toeplitz_prod.hpp
 *   
 *   Fast symmetric Toeplitz matrix-vector product
 *
 *   Copyright (C) 2024   Guilherme F. Fornel        <gffrnl@gmail.com>
 *                        Fabio Souto de Azevedo     <fazedo@gmail.com>
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

#include <cstddef>
#include <cassert>
#include <cmath>
#include <algorithm>
#ifdef FSYMMTP_MANAGER_LOG
    #include <iostream>
#endif
#include <b118/utility.hpp>
#include <b118/numbers.hpp>
#include <b118/fftw.hpp>

namespace b118 { namespace linalg { namespace toeplitz {

    template<typename Real>
    struct fast_symm_toeplitz_prod final {
        template<
            class InputBidirIt1, class InputBidirIt2,
            class OutputBidirIt
        >
        OutputBidirIt operator()(InputBidirIt1 row_beg, InputBidirIt1 row_end,
                                 InputBidirIt2 vec_beg,
                                 OutputBidirIt prod_beg,
                                 bool inplace = false) {
            auto const len = std::distance(row_beg, row_end);

            {
                if (len < static_cast<decltype(len)>(0))
                    throw
                        std::invalid_argument("fast_symm_toeplitz_prod(): "
                                              "std::distance(row_beg, row_end)"
                                              " < 0");

                if (len == static_cast<decltype(len)>(0))
                    return prod_beg;  // TODO(gffrnl): prod_beg  OR
                                      //               std::next(prod_beg)??

                if (len == static_cast<decltype(len)>(1)) {
                    *prod_beg = (*row_beg) * (*vec_beg);
                    return std::next(prod_beg);
                }

                manager::resize(static_cast<std::size_t>(len));
            }

            std::size_t size     = manager::size();
            std::size_t padding  = manager::padding();
            std::size_t aug_size = manager::aug_size();

            Real * const aug_row = manager::aug_row();
            Real * const aug_vec = manager::aug_vec();
            Real * const aux     = manager::aux();

            auto const & plan = manager::plan();

            // Construct the augmented arrays
            std::copy(row_beg, row_end, aug_row);
            for (std::size_t i = 2; i < size; ++i)
                aug_row[i-2] -= (*std::next(row_beg, i));
            std::copy_n(vec_beg, size, aug_vec + padding + 1);

            // Populate with 0 the rest of
            std::fill_n(aug_row + size, aug_size - size, 0);
            std::fill_n(aug_vec, padding + 1, 0);
            std::fill_n(aug_vec + padding + 1 + size,
                        aug_size - (padding + 1 + size), 0);


            // // Populate with 0 the rest of
            // std::memset(static_cast<void *>(aug_mat_first_row + size),
            //             0,
            //             (aug_size - size) * sizeof(double));
            // std::memset(static_cast<void *>(aug_vec),
            //             0,
            //             (padding + 1) * sizeof(double));
            // std::memset(static_cast<void *>(aug_vec + padding + 1 + size),
            //             0,
            //             (aug_size - (padding + 1 + size)) * sizeof(double));

            // Compute the DST1 of aug_mat_first_row and store in aux
            // fftw_execute_r2r(plan.ptr, aug_row, aux);
            // fftw_r2r_wrapper<Real>::execute(plan, aug_row, aux);
            plan.execute(aug_row, aux);

            // Compute the DST1 of y and store in aug_mat_first_row
            // fftw_execute_r2r(plan.ptr, aug_vec, aug_row);
            // fftw_r2r_wrapper<Real>::execute(plan, aug_vec, aug_row);
            plan.execute(aug_vec, aug_row);

            // Multiply
            {
                Real const scaling = static_cast<Real>(1)
                    / static_cast<Real>(4 * (aug_size + 1));
                Real const delta_theta = b118::numbers::pi_v<Real>
                    / static_cast<Real>(aug_size + 1);
                for (std::size_t i = 0; i < aug_size; ++i)
                    aux[i] *= (scaling * aug_row[i]
                            / sin(static_cast<Real>((i + 1) * delta_theta)));
            }
            // NOTE: TODO(gffrnl) REWRITE THIS NOTE
            // need to multiply lamb by 1/2 rather
            //          sqrt(n+1)/2 because in fftw3 scaling
            //          of DST1 is 2

            // Now compute the DST1 of aux
            // fftw_execute_r2r(plan.ptr, aux, aug_row);
            // fftw_r2r_wrapper<Real>::execute(plan, aux, aug_row);
            plan.execute(aux, aug_row);

            // Copy only necessary values to prod
            if (!inplace) {
                std::copy_n(aug_row + padding + 1, size, prod_beg);
            } else {
                for (std::size_t i = 0; i < size; ++i) {
                    // b[i] += mu[i + k + 1];
                    *(prod_beg + i) += *(aug_row + i +padding + 1);
                }
            }

            return std::next(prod_beg, size);
        }

     public:
        class manager final : private b118::fftw::managed<Real> {
            std::size_t capacity_;
            std::size_t size_;
            std::size_t padding_;
            std::size_t aug_size_;

            Real * aug_row_;
            Real * aug_vec_;
            Real * aux_;

            b118::fftw::plan<Real, b118::fftw::r2r> plan_;

         public:
            static manager& initialize(std::size_t n = 0) {
                static manager obj(n);
                return obj;
            }

            static std::size_t capacity() {
                return initialize().capacity_;
            }

            static std::size_t size() {
                return initialize().size_;
            }

            static std::size_t padding() {
                return initialize().padding_;
            }

            static std::size_t aug_size() {
                return initialize().aug_size_;
            }

            static manager & resize(std::size_t size) {
                return initialize().resize_impl(size);
            }

            static Real * const aug_row() {
                return initialize().aug_row_;
            }

            static Real * const aug_vec() {
                return initialize().aug_vec_;
            }

            static Real * const aux() {
                return initialize().aux_;
            }

            static b118::fftw::plan<Real, b118::fftw::r2r> const & plan() {
                return initialize().plan_;
            }

         private:
            // Deleted copy constructor and copy assignment:
            manager(manager &) = delete;
            manager& operator=(manager &) = delete;

            manager(std::size_t n = 2) {
                capacity_ = (n > 2)? n : 2;
                size_     = 0;
                padding_  = 0;
                aug_size_ = 0;
                aug_row_  = nullptr;
                aug_vec_  = nullptr;
                aux_      = nullptr;
                plan_     = b118::fftw::plan<Real, b118::fftw::r2r>();
                if (n >= 2) {
                    set_aug_size(capacity_);
                    alloc_fftw(aug_size_);
                }
            }

            ~manager() {
#ifdef FSYMMTP_MANAGER_LOG
                std::clog << "*** b188::linalg::toeplitz::"
                          << "fast_symm_toeplitz_prod::manager"
                          << " destructor called ***"
                          << std::endl;
#endif
                free_vars();
            }

            void set_aug_size(std::size_t size) {
                padding_  = (size % 2 == 0) ? (size - 4) / 2 : (size - 3) / 2;
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
                //                            fftw_r2r_kind kind,
                //                            unsigned flags);
                //
                // * n gives the (physical) size of the transform dimensions.
                //   They can be any positive integer.
                //
                //   - FFTW is generally best at handling sizes of the form
                //     2^a 3^b 5^c 7^d 11^e 13^f, where e+f is either 0 or 1,
                //     and the other exponents are arbitrary. Other sizes are
                //     computed by means of a slow, general-purpose algorithm 
                //     (which nevertheless retains O(n log n) performance even
                //     for prime sizes). (It is possible to customize FFTW for
                //     different array sizes; see Installation and
                //     Customization.)
                //     Transforms whose sizes are powers of 2 are especially
                //     fast.
                //
                //   - For a REDFT00 or RODFT00 transform kind in a dimension
                //     of size n, it is n-1 or n+1, respectively, that should
                //     be factorizable in the above form.
                //
                // ...
                //
                // (gffrnl) So, since we use RODFT00, we choose the augmented
                // (gffrnl) size of the form 2^a - 1:
                aug_size_ = b118::next_exp2(size + 2 * (padding_ + 1)) - 1;
            }

            manager & resize_impl(std::size_t size) {
                if (size != size_) {
                    std::size_t const old_aug_size = aug_size_;
                    set_aug_size(size);
                    if (aug_size_ > old_aug_size) {
                        // // Remark: first, nullify plan_
                        // plan_ = b118::fftw::plan<Real, b118::fftw::r2r>();
                        // free_vars();
                        alloc_fftw(aug_size_);
                    }
                    // resize the plan
                    resize_plan();
                    size_     = size;
                    capacity_ = size;
                }
                return *this;
            }

            void free_vars() {
                if (aux_ != nullptr) {
                    b118::fftw::free<Real>::real(aux_);
                    aux_ = nullptr;
                }
                if (aug_vec_ != nullptr) {
                    b118::fftw::free<Real>::real(aug_vec_);
                    aug_vec_ = nullptr;
                }
                if (aug_row_ != nullptr) {
                    b118::fftw::free<Real>::real(aug_row_);
                    aug_row_ = nullptr;
                }
            }

            void alloc_fftw(std::size_t n) {
                // Remark: first, nullify plan_
                plan_ = b118::fftw::plan<Real, b118::fftw::r2r>();
                
                // Then, free pointers
                free_vars();
                
                assert(aug_row_ == nullptr);
                assert(aug_vec_ == nullptr);
                assert(aux_     == nullptr);

                aux_     = b118::fftw::alloc<Real>::real(n);
                assert(aux_     != nullptr);

                aug_vec_ = b118::fftw::alloc<Real>::real(n);
                assert(aug_vec_ != nullptr);

                aug_row_ = b118::fftw::alloc<Real>::real(n);
                assert(aug_row_ != nullptr);
            }

            void resize_plan() {
                plan_ = b118::fftw::plan<
                            Real,
                            b118::fftw::r2r
                        >(aug_size_,aug_row_, aug_vec_,
                          //FFTW_RODFT00,
                          b118::fftw::r2rkind::rodft00,
                          FFTW_ESTIMATE);
            }
        };
    };

}}}  // end namespace b118::linalg::toeplitz
