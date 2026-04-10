// Benchmark for proximity solver kernel
// Compile: g++ -std=c++17 -O2 -I cpp/include -o bench_proximity bench_proximity.cpp <sources>

#include "supermag/proximity.h"
#include <chrono>
#include <cstdio>
#include <vector>

int main() {
    const int n_points = 10000;
    const int n_iterations = 1000;
    std::vector<double> x(n_points), F(n_points);

    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < n_iterations; ++iter) {
        supermag_proximity_pair_amplitude(1.0, 2.0, 50.0, n_points, x.data(), F.data());
    }
    auto end = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::printf("pair_amplitude: %d iterations x %d points\n", n_iterations, n_points);
    std::printf("  Total: %.2f ms\n", ms);
    std::printf("  Per call: %.3f us\n", ms * 1000.0 / n_iterations);
    return 0;
}
