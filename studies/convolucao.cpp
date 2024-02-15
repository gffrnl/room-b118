// g++ convolucao.cpp -lfftw3 -lopenblas -lpthread -lbenchmark -O3
#include <fftw3.h>
#include <x86_64-linux-gnu/cblas.h>
// #include <benchmark/benchmark.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <typeinfo>
#include <cassert>
// #include "./fft.h" //  // #define BENCH


std::vector<double> conv1d_ingenua(size_t n0, size_t na, int ja, std::vector<double>& mu, std::vector<double>& y) {
                            
    //for (std::size_t i = 0; i < n0; ++i)
    //    for (std::ptrdiff_t j = 0; j < na; ++j)
    //        Ba[i*na + j] = mu[ja+i-j];

    // frLap_y[i] = produto interno da linha i de Ba com y
    // a linha i de Ba é dada por mu[ja+i-j] com j variando de 0 a na-1

    std::vector<double> frLap_y(n0);
    for (std::size_t i = 0; i < n0; ++i) {
        frLap_y[i] = cblas_ddot(
                        na,
                        mu.data() + ja + i + (-na + 1),
                        -1,      // incremento em mu = -1
                        y.data(),
                        1);
    }
    // Aqui tem um gotcha: se o incremento é negativo, a soma começa do final:
    /*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
    */

    return frLap_y;
}

std::vector<double> cross_ingenua(size_t n0, size_t nb, int jb, std::vector<double>& mu, std::vector<double>& y) {
    
  
    //  for (std::ptrdiff_t i = 0; i < n0; ++i)
    //      for (std::size_t j = 0; j < nb; ++j)
    //          Bb.at(i*nb + j) = mu.at(n0-i+j);

    // frLap_y[i] = produto interno da linha i de Bb com y descolado de jb+1
    // a linha i de Bb é dada por mu[n0-i+j] com j variando de 0 a nb-1
    // frLap_y[i] = cblas_ddot(nb, mu.data() + n0 - i, 1, y.data() + jb + 1, 1);

    std::vector<double> frLap_y(n0);
    for (std::size_t i = 0; i < n0; ++i) {
        frLap_y[i] = cblas_ddot(
                        nb,
                        mu.data() + n0 - i - 1,  // + n0 - i,
                        1,      // incremento em mu = 1
                        y.data(),  // + jb + 1,
                        1);
    }
    return frLap_y;
}
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
                size_t n0,  // TODO(Guilherme/Fabio): nomes mais descritivos
                size_t na,  // size of y. redundant. Deixo isso?
                size_t ja,
                const std::vector<double>& mu,
                const std::vector<double>& y,
                std::vector<double> *output,
                bool cross_correlation = false) {

        // Copy y to in
        if (cross_correlation == true) {
            std::fill(in, in + conv_size - y.size(), 0.0);
            std::copy(y.begin(), y.end(), in + conv_size - y.size());
        } else {
            std::copy(y.begin(), y.end(), in);
            std::fill(in + y.size(), in + conv_size, 0.0);
        }


        // mu to kernel. Note that we are copying na + n0 - 1 elements
        // Hence the kernel is always larger than y, which is not a problem
        // at all, but it may be a bit confusing.
        // Warning: we need to copy only the slice of mu that we are going to use
        // The indices must be calculated before calling this function
        auto ctn = std::copy(mu.begin() + ja - na + 1,
                             mu.begin() + ja + n0,
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

        // Normalize
        // we may reescale only the slice of interest
        // For debugging purposes, we may want to see the full output
        for (size_t i = 0; i < conv_size; ++i) {
            out[i] /= conv_size;
        }

        // Copy output
        output->resize(n0);
        // double* beginning_valid_window = out + std::min(y.size(), n0 + na) - 1;
        double* beginning_valid_window =
                        out + (cross_correlation ? 0 : y.size() - 1);

        std::copy(beginning_valid_window, beginning_valid_window + n0, output->begin());
        return *this;
    }

    void get_full_output(std::vector<double> *v) const {
        if (v->size() != conv_size) {
            v->resize(conv_size);
        }

        std::copy(out, out + conv_size, v->begin());
    }

};

#ifndef BENCH

int main() {
    
    // convolution conv(na + n0);
    double maior_erro_max = 0;
    convolution conv(1 << 16);
    for (int na = 5; na < 16; ++na) {
        for (int n0 = 3; n0 < 14; n0 += 1) {
    
            // int na = 12;  // dimensão da entrada y
            // int n0 = 12;  // dimensão da saída
            int ja = (na-1) + 0;  // índice do primeiro elemento de mu, >= na-1.
            
            // n0 = dimensao da saída
            // na = dimensao da entrada y 
            // Assim temos:
            // frLap_y[i] = sum(mu[ja + i - j] * y[j], j = 0 ... na-1), i = 0 .. n0 - 1
            // range de mu: [ja - na + 1, ja + n0 - 1] -> na + n0

            // A convolução via FFT exige na + (na + n0) - 1 elementos
            // isto é 2*na + n0 - 1
            // O algoritmo ingênuo exige n0 * na operações

            std::vector<double> mu(na + n0 - 1);  // "kernel"
            std::vector<double> y(na);        // entrada
            std::vector<double> frLap_y(n0);  // saída
            

                
            for (int j = 0; j < na; ++j) {
                y.at(j) = j + 1;
            }

            for (int j = 0; j < mu.size(); ++j) {
                mu.at(j) = j*1.0/(j+1);
            }

            
            conv.create_plans(na + n0 - 1);
            
            conv.conv(n0, na, ja, mu, y, &frLap_y, true);
            // conv.get_output(&frLap_y);

            for (int i = 0; i < frLap_y.size(); ++i) {
                std::cout << i << "   :" << std::setprecision(10) << frLap_y.at(i) << '\n';
         }
            std::cout << "----------\n";

            std::vector<double> full_output;
            conv.get_full_output(&full_output);
            for (int i = 0; i < full_output.size(); ++i) {
                std::cout << i << "   :" << std::setprecision(10) << full_output.at(i) << '\n';
            }
            std::cout << "----------\n";


            std::vector<double> frLap_y_2 = cross_ingenua(n0, na, ja,  mu, y);
            // std::vector<double> frLap_y_2 = conv1d_ingenua(n0, na, ja,  mu, y);
            for (int i = 0; i < frLap_y_2.size(); ++i) {
                std::cout << i << "   :" << std::setprecision(10) << frLap_y_2.at(i) << '\n';
         }

            double erro_max = 0;
            for (int i = 0; i < n0; ++i) {
                double fr_1 = frLap_y.at(i);
                double fr_2 = frLap_y_2.at(i);
                double erro = std::fabs(fr_1 - fr_2);  // (std::fabs(fr_1) + std::fabs(fr_2) + 1e-10);
                erro_max = std::max(erro_max, erro);
                if (erro > 1e-10) {
                    std::cout << "na = " << na << " n0 = " << n0 << '\n';
                    std::cout << fr_1 << " : " << fr_2 << " : " << erro << '\n';
                    assert(false);
                }
            }
            maior_erro_max = std::max(maior_erro_max, erro_max);
            std::cout << "erro = " << erro_max << "  maior_erro = " << maior_erro_max << '\n';
        }
    }

    
    return 0;
}

#endif
#ifdef BENCH

static void ingenua(benchmark::State& state) {
    auto N = state.range(0);

    int na = N;  // dimensão da entrada y
    int n0 = N/10;  // dimensão da saída
    int ja = (na-1) + 0;  // índice do primeiro elemento de mu, >= na-1.

    std::vector<double> mu(na + n0 - 1);  // "kernel"
    std::vector<double> y(na);        // entrada
    std::vector<double> frLap_y(n0);  // saída

    for (int j = 0; j < na; ++j) {
        y.at(j) = 2.0/(10.0 + j)*j  + 1;
    }

    for (int j = 0; j < na + n0 - 1; ++j) {
        mu.at(j) = 1.0/(j+1);
    }
    // mu.at(0)= 100;


    for (auto _ : state) {
        frLap_y = conv1d_ingenua(n0, na, ja, mu, y);
        benchmark::DoNotOptimize(frLap_y);  // Prevent compiler optimizations
    }

    std::stringstream ss;
    ss << "y[3] = " << std::setprecision(5) << frLap_y[3];
    state.SetLabel(ss.str());

    state.SetComplexityN(state.range(0));

}

static void conv_fft(benchmark::State& state) {
    auto N = state.range(0);

    int na = N;  // dimensão da entrada y
    int n0 = N/10;  // dimensão da saída
    int ja = (na-1) + 0;  // índice do primeiro elemento de mu, >= na-1.

    std::vector<double> mu(na + n0 - 1);  // "kernel"
    std::vector<double> y(na);        // entrada
    std::vector<double> frLap_y(n0);  // saída

    for (int j = 0; j < na; ++j) {
        y.at(j) = 2.0/(10.0 + j)*j  + 1;
    }

    for (int j = 0; j < na + n0 - 1; ++j) {
        mu.at(j) = 1.0/(j+1);
    }
    // mu.at(na-1)= 100;

    for (auto _ : state) {
        // auto M = 1 << static_cast<int>(std::log2(mu.size() + 2*y.size()-1) + 1);
        auto M = mu.size();  // + y.size()-2;
        convolution conv(M);
        conv.create_plans(M);
        conv.conv(n0, na, ja, mu, y, &frLap_y);
        assert(std::isnan(frLap_y[2]) == false);
        benchmark::DoNotOptimize(frLap_y);  // Prevent compiler optimizations
    }

    std::stringstream ss;
    ss << "y[3] = " << std::setprecision(5) << frLap_y[3];
    state.SetLabel(ss.str());

    state.SetComplexityN(state.range(0));
}


static void conv_fft_2(benchmark::State& state) {
    auto N = state.range(0);

    int na = N;  // dimensão da entrada y
    int n0 = N/10;  // dimensão da saída
    int ja = (na-1) + 0;  // índice do primeiro elemento de mu, >= na-1.

    std::vector<double> mu(na + n0 - 1);  // "kernel"
    std::vector<double> y(na);        // entrada
    std::vector<double> frLap_y(n0);  // saída

    for (int j = 0; j < na; ++j) {
        y.at(j) = 2.0/(10.0 + j)*j  + 1;
    }

    for (int j = 0; j < na + n0 -1 ; ++j) {
        mu.at(j) = 1.0/(j+1);
    }
    // mu.at(na-1)= 100;
    auto M = 1 << static_cast<int>(std::log2(mu.size()) + 1);
    // auto M = mu.size();
    convolution conv(M);
    conv.create_plans(M);
    for (auto _ : state) {
        conv.conv(n0, na, ja, mu, y, &frLap_y);
        assert(std::isnan(frLap_y[2]) == false);
        benchmark::DoNotOptimize(frLap_y);  // Prevent compiler optimizations
    }

    std::stringstream ss;
    ss << "y[3] = " << std::setprecision(5) << frLap_y[3];
    state.SetLabel(ss.str());

    state.SetComplexityN(state.range(0));
}


BENCHMARK(conv_fft)->RangeMultiplier(2)->Range(1<<6, 1<<14)->Complexity(benchmark::oNLogN);  // benchmark::oN);
BENCHMARK(conv_fft_2)->RangeMultiplier(2)->Range(1<<6, 1<<14)->Complexity(benchmark::oNLogN);  // benchmark::oN);
BENCHMARK(ingenua)->RangeMultiplier(2)->Range(1<<6, 1<<14)->Complexity(benchmark::oNSquared);  // benchmark::oN);

BENCHMARK_MAIN();

#endif