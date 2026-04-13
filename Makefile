# SUPERMag — Top-level Makefile
# Builds C++ library, runs tests, builds Python wheels.
#
# Targets:
#   build    — Compile C++ into static library
#   test     — Run C++ unit tests
#   pytest   — Run Python test suite
#   validate — Run validation known-answer tests
#   bench    — Run benchmarks
#   wheel    — Build Python wheel
#   docs     — (placeholder)
#   clean    — Remove all build artifacts

include config/detect_platform.mk

# ── Sources ──────────────────────────────────────────────────
LIB_SRCS = \
	cpp/src/common/error.cpp \
	cpp/src/common/constants.cpp \
	cpp/src/common/digamma.cpp \
	cpp/src/common/allocator.cpp \
	cpp/src/proximity/pair_amplitude.cpp \
	cpp/src/proximity/critical_temp.cpp \
	cpp/src/proximity/kernels.cpp \
	cpp/src/proximity/depairing.cpp \
	cpp/src/proximity/transfer_matrix.cpp \
	cpp/src/proximity/kernel_snf.cpp \
	cpp/src/proximity/kernel_graded.cpp \
	cpp/src/proximity/kernel_domains.cpp \
	cpp/src/linalg/tridiag.cpp \
	cpp/src/linalg/simd_kernels.cpp \
	cpp/src/linalg/eigen.cpp \
	cpp/src/solvers/usadel.cpp \
	cpp/src/solvers/eilenberger.cpp \
	cpp/src/solvers/bdg.cpp \
	cpp/src/solvers/ginzburg_landau.cpp \
	cpp/src/solvers/josephson.cpp \
	cpp/src/solvers/triplet.cpp \
	cpp/src/solvers/root_scalar.cpp \
	cpp/src/solvers/determinant.cpp \
	cpp/src/proximity/spin_active.cpp \
	cpp/src/proximity/depairing_models.cpp \
	cpp/src/proximity/optimizer.cpp \
	cpp/src/proximity/interface.cpp

LIB_OBJS = $(patsubst cpp/%.cpp,build/obj/%$(OBJ_EXT),$(LIB_SRCS))
STATIC_LIB = build/$(LIB_PREFIX)supermag$(LIB_EXT)

TEST_SRCS = cpp/test/test_proximity.cpp cpp/test/test_tridiag.cpp cpp/test/test_stubs.cpp \
	cpp/test/test_digamma.cpp cpp/test/test_root_scalar.cpp cpp/test/test_determinant.cpp \
	cpp/test/test_transfer_matrix.cpp cpp/test/test_kernel_snf.cpp \
	cpp/test/test_kernel_graded.cpp cpp/test/test_kernel_domains.cpp \
	cpp/test/test_usadel.cpp cpp/test/test_eilenberger.cpp \
	cpp/test/test_bdg.cpp cpp/test/test_ginzburg_landau.cpp \
	cpp/test/test_josephson.cpp cpp/test/test_triplet.cpp \
	cpp/test/test_depairing_models.cpp cpp/test/test_optimizer.cpp
TEST_BINS = $(patsubst cpp/test/%.cpp,build/test/%$(EXE_EXT),$(TEST_SRCS))

# ── Build ────────────────────────────────────────────────────
.PHONY: build test pytest validate bench wheel docs clean shared ocaml ocaml-test

build: $(STATIC_LIB)

ifeq ($(OS),Windows_NT)
MKDIR_P = if not exist "$(subst /,\,$1)" mkdir "$(subst /,\,$1)"
else
MKDIR_P = mkdir -p $1
endif

build/obj/%$(OBJ_EXT): cpp/%.cpp
	@$(call MKDIR_P,$(dir $@))
	$(call compile_obj)

$(STATIC_LIB): $(LIB_OBJS)
	@$(call MKDIR_P,$(dir $@))
	$(call archive_lib)

# ── Shared Library (for ctypes FFI) ──────────────────────────
ifeq ($(OS),Windows_NT)
SHARED_LIB = build/supermag.dll
$(SHARED_LIB): $(LIB_OBJS)
	@$(call MKDIR_P,$(dir $@))
	link /DLL /OUT:$@ $^
else
SHARED_LIB = build/libsupermag.so
$(SHARED_LIB): $(LIB_OBJS)
	@$(call MKDIR_P,$(dir $@))
	$(CXX) -shared -o $@ $^ -lm
endif

shared: $(SHARED_LIB)

# ── C++ Tests ────────────────────────────────────────────────
test: $(TEST_BINS)
	@echo "Running C++ tests..."
ifeq ($(OS),Windows_NT)
	$(foreach t,$(TEST_BINS),$(subst /,\,$(t)) &&) echo "All C++ tests passed."
else
	@for t in $(TEST_BINS); do \
		echo "--- $$t ---"; \
		$$t || exit 1; \
	done
	@echo "All C++ tests passed."
endif

build/test/%$(EXE_EXT): cpp/test/%.cpp $(STATIC_LIB)
	@$(call MKDIR_P,$(dir $@))
	$(call link_exe)

# ── Python ───────────────────────────────────────────────────
pytest:
	python -m pytest python/tests/ -v

validate:
	python validation/run_validation.py

bench:
	python benchmarks/bench_runner.py

wheel:
	pip wheel python/ -w dist/

docs:
	@echo "Documentation generation not yet configured."
	@echo "See docs/ for existing theory and tutorial documents."

# ── OCaml ────────────────────────────────────────────────────
ocaml: shared
	cd ocaml && dune build

ocaml-test: ocaml
	cd ocaml && dune runtest

# ── Clean ────────────────────────────────────────────────────
ifeq ($(OS),Windows_NT)
clean:
	if exist build rmdir /s /q build
	if exist dist rmdir /s /q dist
	if exist _skbuild rmdir /s /q _skbuild
	for /r . %%d in (__pycache__) do @if exist "%%d" rmdir /s /q "%%d"
	del /s /q *.pyc 2>nul || ver >nul
else
clean:
	rm -rf build/ dist/ _skbuild/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
endif
