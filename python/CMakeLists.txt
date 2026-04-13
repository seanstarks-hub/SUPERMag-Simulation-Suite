cmake_minimum_required(VERSION 3.15)
project(supermag_python)

# Find required packages
find_package(pybind11 CONFIG REQUIRED)

# Include directories for SUPERMag headers
include_directories(${CMAKE_SOURCE_DIR}/../cpp/include)
include_directories(${CMAKE_SOURCE_DIR}/../cpp/src)

# Find the C library
# Look for the pre-built static library
find_library(SUPERMAG_LIBRARY 
    NAMES supermag
    PATHS ${CMAKE_SOURCE_DIR}/../build
          ${CMAKE_BINARY_DIR}/../build
    NO_DEFAULT_PATH
    REQUIRED
)

# Create the Python module
pybind11_add_module(_native supermag/_binding.cpp)

# Link against the SUPERMag C library
target_link_libraries(_native PRIVATE ${SUPERMAG_LIBRARY})

# Install the extension into the supermag package
install(TARGETS _native DESTINATION supermag)
