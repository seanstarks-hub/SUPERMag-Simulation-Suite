# SUPERMag — Top-level Makefile
# Dev-convenience wrapper around CMake.  The canonical build entry point is
# python/CMakeLists.txt; this Makefile simply invokes it with useful presets.
#
# Targets:
#   build    — Compile C++ static library via CMake
#   shared   — Build shared library for OCaml ctypes FFI
#   test     — Compile & run C++ unit tests
#   pytest   — Run Python test suite
#   validate — Run validation known-answer tests
#   bench    — Run benchmarks
#   wheel    — Build Python wheel
#   docs     — (placeholder)
#   ocaml    — Build OCaml project (depends on shared)
#   ocaml-test — Run OCaml tests
#   clean    — Remove all build artifacts

CXX ?= g++

.PHONY: build shared test pytest validate bench wheel docs ocaml ocaml-test clean

# ── Build ────────────────────────────────────────────────────
build:
	cmake -B build -S python -DCMAKE_BUILD_TYPE=Release
	cmake --build build

shared:
	cmake -B build -S python -DCMAKE_BUILD_TYPE=Release -DSUPERMAG_BUILD_SHARED=ON
	cmake --build build
	cmake --build build

# ── C++ Tests ────────────────────────────────────────────────
test: build
	@echo "Running C++ tests..."
	@mkdir -p build/test
	@for src in cpp/test/test_*.cpp; do \
		name=$$(basename $$src .cpp); \
		$(CXX) -std=c++17 -O2 -I cpp/include -I cpp/src -o build/test/$$name $$src -Lbuild -lsupermag_core -lm && \
		echo "--- $$name ---" && \
		build/test/$$name || exit 1; \
	done
	@echo "All C++ tests passed."

# ── Python ───────────────────────────────────────────────────
pytest:
	python -m pytest python/tests/ -v

validate:
	python validation/run_validation.py

bench:
	python benchmarks/bench_runner.py

wheel:
	pip wheel . -w dist/

docs:
	@echo "Documentation generation not yet configured."
	@echo "See docs/ for existing theory and tutorial documents."

# ── OCaml ────────────────────────────────────────────────────
ocaml: shared
	cd ocaml && dune build

ocaml-test: ocaml
	cd ocaml && dune runtest

# ── Clean ────────────────────────────────────────────────────
clean:
	rm -rf build/ dist/ _skbuild/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
