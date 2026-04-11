// AVX2 intrinsic kernels for hot inner loops
// These provide SIMD-accelerated versions of common operations.
//
// Compile with: /arch:AVX2 (MSVC) or -mavx2 (GCC/Clang)
// At runtime, only called if AVX2 is detected.

#include <cstddef>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

extern "C" {

/**
 * SIMD-accelerated: out[i] = a[i] * b[i] for i in [0, n)
 * Falls back to scalar if AVX2 not available at compile time.
 */
void supermag_simd_vec_mul(const double* a, const double* b, double* out, int n) {
    int i = 0;
#if defined(__AVX2__)
    for (; i + 3 < n; i += 4) {
        __m256d va = _mm256_loadu_pd(a + i);
        __m256d vb = _mm256_loadu_pd(b + i);
        __m256d vr = _mm256_mul_pd(va, vb);
        _mm256_storeu_pd(out + i, vr);
    }
#endif
    for (; i < n; ++i) {
        out[i] = a[i] * b[i];
    }
}

/**
 * SIMD-accelerated: out[i] += scale * x[i] (DAXPY-like)
 */
void supermag_simd_axpy(double scale, const double* x, double* out, int n) {
    int i = 0;
#if defined(__AVX2__)
    __m256d vs = _mm256_set1_pd(scale);
    for (; i + 3 < n; i += 4) {
        __m256d vx = _mm256_loadu_pd(x + i);
        __m256d vo = _mm256_loadu_pd(out + i);
        vo = _mm256_add_pd(vo, _mm256_mul_pd(vs, vx));
        _mm256_storeu_pd(out + i, vo);
    }
#endif
    for (; i < n; ++i) {
        out[i] += scale * x[i];
    }
}

}
