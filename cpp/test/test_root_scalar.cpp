// Unit tests for the root_scalar solver (Brent's method with grid scan)
#define _USE_MATH_DEFINES
#include "../src/solvers/root_scalar.h"
#include <cstdio>
#include <cmath>
#include <cassert>

// f(x) = x^2 - 4  → roots at x = -2 and x = 2
static double quadratic(double x, void*) {
    return x * x - 4.0;
}

// f(x) = sin(x) → roots at 0, pi, 2pi, ...
static double sine_func(double x, void*) {
    return std::sin(x);
}

// f(x) = (x - 1)(x - 3)(x - 7)  → roots at 1, 3, 7
static double cubic(double x, void*) {
    return (x - 1.0) * (x - 3.0) * (x - 7.0);
}

// f(x) = x - 5  → single root at x = 5
static double linear(double x, void*) {
    return x - 5.0;
}

void test_quadratic_highest() {
    // Should return 2.0 (highest root in [0, 10])
    double root = supermag::root_scalar_solve(quadratic, nullptr, 0.0, 10.0, 1e-12);
    assert(std::abs(root - 2.0) < 1e-9);
    std::printf("  PASS: test_quadratic_highest\n");
}

void test_quadratic_both() {
    // In [-10, 10], roots at -2 and 2. Highest = 2.
    double root = supermag::root_scalar_solve(quadratic, nullptr, -10.0, 10.0, 1e-12);
    assert(std::abs(root - 2.0) < 1e-9);
    std::printf("  PASS: test_quadratic_both\n");
}

void test_sine_highest() {
    // sin(x) in [0.1, 10]: roots at pi, 2pi, 3pi. Highest = 3pi ≈ 9.42
    double root = supermag::root_scalar_solve(sine_func, nullptr, 0.1, 10.0, 1e-12);
    assert(std::abs(root - 3.0 * M_PI) < 1e-6);
    std::printf("  PASS: test_sine_highest\n");
}

void test_cubic_highest() {
    // Roots at 1, 3, 7. Highest = 7.
    double root = supermag::root_scalar_solve(cubic, nullptr, 0.0, 10.0, 1e-12);
    assert(std::abs(root - 7.0) < 1e-9);
    std::printf("  PASS: test_cubic_highest\n");
}

void test_linear() {
    double root = supermag::root_scalar_solve(linear, nullptr, 0.0, 10.0, 1e-12);
    assert(std::abs(root - 5.0) < 1e-9);
    std::printf("  PASS: test_linear\n");
}

void test_no_root() {
    // f(x) = x^2 + 1 has no real roots
    auto f = [](double x, void*) -> double { return x * x + 1.0; };
    double root = supermag::root_scalar_solve(f, nullptr, -10.0, 10.0, 1e-12);
    assert(std::isnan(root));
    std::printf("  PASS: test_no_root\n");
}

int main() {
    std::printf("Running root_scalar tests...\n");
    test_quadratic_highest();
    test_quadratic_both();
    test_sine_highest();
    test_cubic_highest();
    test_linear();
    test_no_root();
    std::printf("All root_scalar tests passed!\n");
    return 0;
}
