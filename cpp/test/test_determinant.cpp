// Unit tests for the determinant root solver
// (Same algorithm as root_scalar, applied to determinant functions)
#define _USE_MATH_DEFINES
#include "../src/solvers/root_scalar.h"
#include <cstdio>
#include <cmath>
#include <cassert>

namespace supermag {
double determinant_solve(double (*f)(double, void*), void*, double, double, double);
}

// Simple "determinant" that is zero at x = 3.5
static double det_linear(double x, void*) {
    return x - 3.5;
}

// 2x2 matrix M(T) = [[T-1, 0.5], [0.5, T-4]]
// det M = (T-1)(T-4) - 0.25 = T^2 - 5T + 3.75
// roots: T = (5 ± sqrt(25-15))/2 = (5 ± sqrt(10))/2
// T ≈ 0.919 and T ≈ 4.081
static double det_2x2(double T, void*) {
    return (T - 1.0) * (T - 4.0) - 0.25;
}

// Oscillatory determinant: det = cos(x) - 0.5
// roots where cos(x) = 0.5, i.e. x = pi/3, 5pi/3, 7pi/3, ...
static double det_oscillatory(double x, void*) {
    return std::cos(x) - 0.5;
}

void test_det_linear() {
    double root = supermag::determinant_solve(det_linear, nullptr, 0.0, 10.0, 1e-12);
    assert(std::abs(root - 3.5) < 1e-9);
    std::printf("  PASS: test_det_linear\n");
}

void test_det_2x2_highest() {
    // Highest root ≈ 4.081
    double expected = (5.0 + std::sqrt(10.0)) / 2.0;
    double root = supermag::determinant_solve(det_2x2, nullptr, 0.0, 10.0, 1e-12);
    assert(std::abs(root - expected) < 1e-6);
    std::printf("  PASS: test_det_2x2_highest\n");
}

void test_det_oscillatory() {
    // In [0, 8], cos(x) = 0.5 at pi/3 ≈ 1.047, 5pi/3 ≈ 5.236, 7pi/3 ≈ 7.330
    // Highest root ≈ 7.330
    double root = supermag::determinant_solve(det_oscillatory, nullptr, 0.0, 8.0, 1e-12);
    assert(std::abs(root - 7.0 * M_PI / 3.0) < 1e-4);
    std::printf("  PASS: test_det_oscillatory\n");
}

int main() {
    std::printf("Running determinant tests...\n");
    test_det_linear();
    test_det_2x2_highest();
    test_det_oscillatory();
    std::printf("All determinant tests passed!\n");
    return 0;
}
