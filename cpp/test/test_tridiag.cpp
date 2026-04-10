// Unit tests for tridiagonal solver
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

extern "C" int supermag_tridiag_solve(
    const double* a, const double* b, const double* c, double* d, int n);

void test_tridiag_3x3() {
    // System: [2 -1 0; -1 2 -1; 0 -1 2] x = [1; 0; 1]
    // Solution: x = [1; 1; 1]
    double a[] = {0, -1, -1};
    double b[] = {2,  2,  2};
    double c[] = {-1, -1, 0};
    double d[] = {1, 0, 1};
    int rc = supermag_tridiag_solve(a, b, c, d, 3);
    assert(rc == SUPERMAG_OK);
    assert(std::abs(d[0] - 1.0) < 1e-10);
    assert(std::abs(d[1] - 1.0) < 1e-10);
    assert(std::abs(d[2] - 1.0) < 1e-10);
    std::printf("  PASS: test_tridiag_3x3\n");
}

void test_tridiag_1x1() {
    double a[] = {0};
    double b[] = {3.0};
    double c[] = {0};
    double d[] = {6.0};
    int rc = supermag_tridiag_solve(a, b, c, d, 1);
    assert(rc == SUPERMAG_OK);
    assert(std::abs(d[0] - 2.0) < 1e-10);
    std::printf("  PASS: test_tridiag_1x1\n");
}

void test_tridiag_null() {
    int rc = supermag_tridiag_solve(nullptr, nullptr, nullptr, nullptr, 3);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_tridiag_null\n");
}

int main() {
    std::printf("Running tridiag tests...\n");
    test_tridiag_3x3();
    test_tridiag_1x1();
    test_tridiag_null();
    std::printf("All tridiag tests passed!\n");
    return 0;
}
