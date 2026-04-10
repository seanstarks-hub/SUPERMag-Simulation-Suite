# GCC/Clang flags for Linux x86_64
CXX      = g++
CXXFLAGS = -std=c++17 -O2 -mavx2 -Wall -Wextra -fPIC
CPPFLAGS = -I cpp/include
AR       = ar
ARFLAGS  = rcs
OBJ_EXT  = .o
LIB_EXT  = .a
LIB_PREFIX = lib
EXE_EXT  =

# Compile rule: source -> object
define compile_obj
$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<
endef

# Archive rule: objects -> static library
define archive_lib
$(AR) $(ARFLAGS) $@ $^
endef

# Link rule: objects -> executable
define link_exe
$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -lm
endef
