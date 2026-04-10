# MSVC compiler flags for Windows x86_64
CXX      = cl
CXXFLAGS = /std:c++17 /O2 /arch:AVX2 /EHsc /W4
CPPFLAGS = /I cpp/include
AR       = lib
ARFLAGS  = /OUT:
OBJ_EXT  = .obj
LIB_EXT  = .lib
LIB_PREFIX =
EXE_EXT  = .exe

# Compile rule: source -> object
define compile_obj
$(CXX) $(CXXFLAGS) $(CPPFLAGS) /c /Fo$@ $<
endef

# Archive rule: objects -> static library
define archive_lib
$(AR) $(ARFLAGS)$@ $^
endef

# Link rule: objects -> executable
define link_exe
$(CXX) $(CXXFLAGS) $(CPPFLAGS) /Fe$@ $^
endef
