// g++ convolucao.cpp -lfftw3 -lopenblas -lpthread -lbenchmark -O3
#include <b118/frlap/gdm/convolucao.hpp>
#include <x86_64-linux-gnu/cblas.h>
#include <benchmark/benchmark.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <typeinfo>
#include <cassert>
// #include "./fft.h" //  // #define BENCH

int main() {


    return 0;
}



// static void conv_fft(benchmark::State& state) {
//     auto N = state.range(0);


//     for (auto _ : state) {
// //
//         benchmark::DoNotOptimize(frLap_y);  // Prevent compiler optimizations
//     }

//     std::stringstream ss;
//     ss << "y[3] = " << std::setprecision(5) << frLap_y[3];
//     state.SetLabel(ss.str());

//     state.SetComplexityN(state.range(0));
// }



// BENCHMARK(conv_fft)->RangeMultiplier(2)->Range(1<<6, 1<<14)->Complexity(benchmark::oNLogN);  // benchmark::oN);
// BENCHMARK(conv_fft_2)->RangeMultiplier(2)->Range(1<<6, 1<<14)->Complexity(benchmark::oNLogN);  // benchmark::oN);
// BENCHMARK(ingenua)->RangeMultiplier(2)->Range(1<<6, 1<<14)->Complexity(benchmark::oNSquared);  // benchmark::oN);

// BENCHMARK_MAIN();

// #endif