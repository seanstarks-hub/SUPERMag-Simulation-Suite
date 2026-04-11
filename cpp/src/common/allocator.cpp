/**
 * SIMD-aligned memory allocator for hot buffers.
 *
 * Provides 32-byte aligned allocation (AVX2 requirement) with RAII-safe
 * free. Falls back to standard malloc on platforms without aligned_alloc.
 */

#include <cstdlib>
#include <cstddef>

#ifdef _MSC_VER
#include <malloc.h>
#endif

extern "C" {

void* supermag_aligned_alloc(size_t alignment, size_t size) {
    if (size == 0) return nullptr;
#ifdef _MSC_VER
    return _aligned_malloc(size, alignment);
#else
    return aligned_alloc(alignment, size);
#endif
}

void supermag_aligned_free(void* ptr) {
#ifdef _MSC_VER
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

} // extern "C"
