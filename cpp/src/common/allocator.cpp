// STUB — Aligned memory allocator for SIMD
// TODO: Implement aligned_alloc / _aligned_malloc wrappers for AVX2-aligned buffers.
//       Must support both MSVC (_aligned_malloc/_aligned_free) and POSIX (aligned_alloc/free).

#include <cstdlib>
#include <cstddef>

namespace supermag {
namespace detail {

// TODO: Implement cross-platform aligned allocation
// void* aligned_alloc(size_t alignment, size_t size);
// void  aligned_free(void* ptr);

} // namespace detail
} // namespace supermag
