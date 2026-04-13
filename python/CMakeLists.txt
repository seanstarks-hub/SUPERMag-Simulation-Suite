cmake_minimum_required(VERSION 3.15)
project(supermag_python LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(SUPERMAG_BUILD_SHARED "Build shared library for ctypes/OCaml FFI" OFF)

# ── Paths ────────────────────────────────────────────────────
set(CPP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../cpp")

# ── Static core library ─────────────────────────────────────
add_library(supermag_core STATIC
    ${CPP_ROOT}/src/common/error.cpp
    ${CPP_ROOT}/src/common/constants.cpp
    ${CPP_ROOT}/src/common/digamma.cpp
    ${CPP_ROOT}/src/common/allocator.cpp
    ${CPP_ROOT}/src/proximity/pair_amplitude.cpp
    ${CPP_ROOT}/src/proximity/critical_temp.cpp
    ${CPP_ROOT}/src/proximity/kernels.cpp
    ${CPP_ROOT}/src/proximity/depairing.cpp
    ${CPP_ROOT}/src/proximity/transfer_matrix.cpp
    ${CPP_ROOT}/src/proximity/kernel_snf.cpp
    ${CPP_ROOT}/src/proximity/kernel_graded.cpp
    ${CPP_ROOT}/src/proximity/kernel_domains.cpp
    ${CPP_ROOT}/src/linalg/tridiag.cpp
    ${CPP_ROOT}/src/linalg/simd_kernels.cpp
    ${CPP_ROOT}/src/linalg/eigen.cpp
    ${CPP_ROOT}/src/solvers/usadel.cpp
    ${CPP_ROOT}/src/solvers/eilenberger.cpp
    ${CPP_ROOT}/src/solvers/bdg.cpp
    ${CPP_ROOT}/src/solvers/ginzburg_landau.cpp
    ${CPP_ROOT}/src/solvers/josephson.cpp
    ${CPP_ROOT}/src/solvers/triplet.cpp
    ${CPP_ROOT}/src/solvers/root_scalar.cpp
    ${CPP_ROOT}/src/solvers/determinant.cpp
    ${CPP_ROOT}/src/proximity/spin_active.cpp
    ${CPP_ROOT}/src/proximity/depairing_models.cpp
    ${CPP_ROOT}/src/proximity/optimizer.cpp
    ${CPP_ROOT}/src/proximity/interface.cpp
)

target_include_directories(supermag_core PUBLIC
    ${CPP_ROOT}/include
    ${CPP_ROOT}/src
)

# ── Optional AVX2 ────────────────────────────────────────────
include(CheckCXXCompilerFlag)
if(MSVC)
    check_cxx_compiler_flag("/arch:AVX2" COMPILER_SUPPORTS_AVX2)
    if(COMPILER_SUPPORTS_AVX2)
        target_compile_options(supermag_core PRIVATE /arch:AVX2)
    endif()
else()
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    if(COMPILER_SUPPORTS_AVX2)
        target_compile_options(supermag_core PRIVATE -mavx2)
    endif()
endif()

# ── Math library (non-Windows) ───────────────────────────────
if(NOT WIN32)
    target_link_libraries(supermag_core PUBLIC m)
endif()

# ── Shared library for ctypes / OCaml FFI ────────────────────
if(SUPERMAG_BUILD_SHARED)
    add_library(supermag_shared SHARED $<TARGET_OBJECTS:supermag_core>)
    set_target_properties(supermag_shared PROPERTIES OUTPUT_NAME supermag)
    target_include_directories(supermag_shared PUBLIC
        ${CPP_ROOT}/include
        ${CPP_ROOT}/src
    )
    if(NOT WIN32)
        target_link_libraries(supermag_shared PUBLIC m)
    endif()
endif()

# ── pybind11 extension (optional — present via scikit-build) ─
find_package(pybind11 CONFIG QUIET)
if(pybind11_FOUND)
    pybind11_add_module(_native supermag/_binding.cpp)
    target_link_libraries(_native PRIVATE supermag_core)
    install(TARGETS _native DESTINATION supermag)
endif()
