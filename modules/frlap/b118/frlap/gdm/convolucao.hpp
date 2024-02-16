#include <fftw3.h>
#include <x86_64-linux-gnu/cblas.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <typeinfo>
#include <cassert>


//
// The input is y and mu, and the output is frLap_y.
// mu must have at least na + n0 - 1 elements
// y must have at least na elements
// The convolution is given by:
// frLap_y[i] = sum(mu[(ja + i - j] * y[j], j = 0 ... na - 1), i = 0 .. n0 - 1
// The cross-correlation is given by:
// frLap_y[i] = sum(mu[n0 - 1 - i + j] * y[j + jb], j = 0 ... nb - 1), i = 0 .. n0 - 1
// Atenção ao colocar no código, é mu[n0 - i + j]
// TODO(Guilherme/Fabio): rever os índices.
class convolution {
    double *in;
    double *kernel;
    double *out;
    fftw_complex* in_fft;
    fftw_complex* kernel_fft;
    fftw_plan plan_r2c = nullptr;
    fftw_plan plan_c2r = nullptr;
    size_t conv_size = 0;
    size_t max_size = 0;

 public:
    convolution(size_t N) : max_size(N) {  // NOLINT
        // N is the logical size, i.e. the number of real samples
        // N must be at least na + n0 - 1
        // We may initialize with a large N and then
        // redo the plans with the actual size of the convolution.
        // This avoids the overhead of allocation and reallocating memory.
        in = static_cast<double*>(fftw_malloc(sizeof(double) * N));
        kernel = static_cast<double*>(fftw_malloc(sizeof(double) * N));
        out = static_cast<double*>(fftw_malloc(sizeof(double) * N));
        in_fft = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * (N/2 + 1)));
        kernel_fft = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * (N/2 + 1)));
    }

    ~convolution() {
        fftw_destroy_plan(plan_r2c);
        fftw_destroy_plan(plan_c2r);
        fftw_free(in);
        fftw_free(kernel);
        fftw_free(out);
        fftw_free(in_fft);
        fftw_free(kernel_fft);
    }

    convolution& create_plans(size_t conv_size) {
        if (plan_r2c != nullptr) {
            fftw_destroy_plan(plan_r2c);
        }
        if (plan_c2r != nullptr) {
            fftw_destroy_plan(plan_c2r);
        }
        // conv_size the number of real samples.
        plan_r2c = fftw_plan_dft_r2c_1d(conv_size, in, in_fft, FFTW_ESTIMATE);
        plan_c2r = fftw_plan_dft_c2r_1d(conv_size, in_fft, out, FFTW_ESTIMATE);
        this->conv_size = conv_size;
        return *this;
    }
    void elementwise_product(bool conj = false) {
        if (conj == false) {
            for (size_t i = 0; i < conv_size/2 + 1; ++i) {
                double re = +in_fft[i][0] * kernel_fft[i][0] - in_fft[i][1] * kernel_fft[i][1];
                double im = +in_fft[i][0] * kernel_fft[i][1] + in_fft[i][1] * kernel_fft[i][0];
                in_fft[i][0] = re;
                in_fft[i][1] = im;
            }
        } else {
            for (size_t i = 0; i < conv_size/2 + 1; ++i) {
                double re = +in_fft[i][0] * kernel_fft[i][0] + in_fft[i][1] * kernel_fft[i][1];
                double im = -in_fft[i][0] * kernel_fft[i][1] + in_fft[i][1] * kernel_fft[i][0];
                in_fft[i][0] = re;
                in_fft[i][1] = im;
            }
        }
    }
    convolution& conv(
                size_t output_size,  // TODO(Guilherme/Fabio): nomes mais descritivos
                size_t input_size,  // size of y. redundant. Deixo isso?
               // size_t ja,
                double const * mu,
                double const * input,
                double *output,
                bool cross_correlation = false) {

        
        // Copy y to in
        if (cross_correlation == true) {
            std::fill(in, in + conv_size - input_size, 0.0);
            std::copy(input, input + input_size, in + conv_size - input_size);
        } else {
            std::copy(input, input + input_size, in);
            std::fill(in + input_size, in + conv_size, 0.0);
        }


        // mu to kernel. Note that we are copying na + n0 - 1 elements
        // Hence the kernel is always larger than y, which is not a problem
        // at all, but it may be a bit confusing.
        // Warning: we need to copy only the slice of mu that we are going to use
        // The indices must be calculated before calling this function
        auto ctn = std::copy(mu,
                             mu + input_size + output_size - 1,
                             kernel);

        std::fill(ctn, kernel + conv_size, 0.0);

        // assert(conv_size <= max_size);
        // Execute plan_r2c for each
        fftw_execute_dft_r2c(plan_r2c, in, in_fft);
        fftw_execute_dft_r2c(plan_r2c, kernel, kernel_fft);

        // Elementwise product
        elementwise_product(cross_correlation);

        // Execute plan_c2r, i.e., inverse transform
        // From doc: the c2r transform destroys its input array even for out-of-place transforms.
        fftw_execute_dft_c2r(plan_c2r, in_fft, out);

        // Copy output
        double* beginning_valid_window =
                        out + (cross_correlation ? 0 : input_size - 1);

        for (size_t i = 0; i < output_size; ++i) {
            output[i] += beginning_valid_window[i]/conv_size;
        }
        return *this;
    }

    // void get_full_output(std::vector<double> *v) const {
    //     if (v->size() != conv_size) {
    //         v->resize(conv_size);
    //     }

    //     std::copy(out, out + conv_size, v->begin());
    // }

};
