VERSION = 0.1.3

CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O3

HEADERS = $(shell find include -name "*.hpp")

BUILD_DIR = build
UNIT_BUILD_DIR = $(BUILD_DIR)/unit
INTEGRATION_BUILD_DIR = $(BUILD_DIR)/integration
BENCH_BUILD_DIR = $(BUILD_DIR)/bench

# High-level targets
# - all: build all unit/integration tests + benchmarks + examples + ms + bi
# - test: build unit/integration tests and run them, then ms and bi app tests
# - bench: build benchmarks
# - examples: build examples
# - clean: clean build files
# - debug: library/tests/example via CXXFLAGS += -g -O0 then all; ms/bi via subdir `make debug`

help:
	@echo "Usage: make [target]"
	@echo "Targets:"
	@echo "  all: build unit/integration tests + benchmarks + examples + ms + bi"
	@echo "  examples: build examples (example target)"
	@echo "  test: run library tests then ./src/ms/ms and ./src/bi/bi tests"
	@echo "  bench: build benchmarks"
	@echo "  clean: clean build files"
	@echo "  debug: debug build (root tests/example + src/ms + src/bi)"

UNIT_TESTS = $(UNIT_BUILD_DIR)/packed_vector_test \
             $(UNIT_BUILD_DIR)/alphabet_test \
             $(UNIT_BUILD_DIR)/columns_test \
             $(UNIT_BUILD_DIR)/api_test \
             $(UNIT_BUILD_DIR)/common_test \
             $(UNIT_BUILD_DIR)/interval_encoding_test \
             $(UNIT_BUILD_DIR)/move_table_test \
             $(UNIT_BUILD_DIR)/move_splitting_test \
             $(UNIT_BUILD_DIR)/move_structure_test \
             $(UNIT_BUILD_DIR)/move_test \
             $(UNIT_BUILD_DIR)/permutation_test \
             $(UNIT_BUILD_DIR)/permutation_random_test \
             $(UNIT_BUILD_DIR)/rlbwt_row_test \
             $(UNIT_BUILD_DIR)/rlbwt_structure_test \
             $(UNIT_BUILD_DIR)/permutation_lf_fl_test \
             $(UNIT_BUILD_DIR)/permutation_phi_phi_inv_test \
             $(UNIT_BUILD_DIR)/rlbwt_interval_encoding_test
INTEGRATION_TESTS = $(INTEGRATION_BUILD_DIR)/rlbwt_test \
                    $(INTEGRATION_BUILD_DIR)/move_structure_test \
                    $(INTEGRATION_BUILD_DIR)/move_test \
                    $(INTEGRATION_BUILD_DIR)/permutation_test
BENCH_TESTS = $(BENCH_BUILD_DIR)/move_bench \
              $(BENCH_BUILD_DIR)/permutation_bench \
              $(BENCH_BUILD_DIR)/rlbwt_bench

.PHONY: ms bi example all test bench clean debug help

all: $(UNIT_TESTS) $(INTEGRATION_TESTS) $(BENCH_TESTS) example ms bi

example: examples/example.cpp examples/example.hpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

ms:
	$(MAKE) -C src/ms

bi:
	$(MAKE) -C src/bi

test: $(UNIT_TESTS) $(INTEGRATION_TESTS) ms bi
	$(UNIT_BUILD_DIR)/packed_vector_test
	$(UNIT_BUILD_DIR)/alphabet_test
	$(UNIT_BUILD_DIR)/columns_test
	$(UNIT_BUILD_DIR)/interval_encoding_test
	$(UNIT_BUILD_DIR)/common_test
	$(UNIT_BUILD_DIR)/api_test
	$(UNIT_BUILD_DIR)/move_table_test
	$(UNIT_BUILD_DIR)/move_splitting_test
	$(UNIT_BUILD_DIR)/move_structure_test
	$(UNIT_BUILD_DIR)/move_test
	$(UNIT_BUILD_DIR)/permutation_test
	$(UNIT_BUILD_DIR)/permutation_random_test
	$(UNIT_BUILD_DIR)/rlbwt_row_test
	$(UNIT_BUILD_DIR)/rlbwt_structure_test
	$(UNIT_BUILD_DIR)/permutation_lf_fl_test
	$(UNIT_BUILD_DIR)/permutation_phi_phi_inv_test
	$(UNIT_BUILD_DIR)/rlbwt_interval_encoding_test
	$(INTEGRATION_BUILD_DIR)/rlbwt_test
	$(INTEGRATION_BUILD_DIR)/move_structure_test
	$(INTEGRATION_BUILD_DIR)/move_test
	$(INTEGRATION_BUILD_DIR)/permutation_test
	./src/ms/ms test data
	./src/bi/bi test data
	@echo "==============================================="
	@echo "     All unit, integration, and app tests passed"
	@echo "==============================================="

bench: $(BENCH_TESTS)

# Unit tests (header-only data structures)
$(UNIT_BUILD_DIR)/packed_vector_test: ./tests/unit/ds/packed_vector_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/alphabet_test: ./tests/unit/ds/alphabet_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/common_test: ./tests/unit/common_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/api_test: ./tests/unit/api_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/columns_test: ./tests/unit/columns_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/interval_encoding_test: ./tests/unit/move/interval_encoding_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_table_test: ./tests/unit/move/move_table_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_splitting_test: ./tests/unit/move/move_splitting_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_structure_test: ./tests/unit/move/move_structure_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_test: ./tests/unit/perm/move_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/permutation_test: ./tests/unit/perm/permutation_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/permutation_random_test: ./tests/unit/perm/permutation_random_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/rlbwt_row_test: ./tests/unit/rlbwt/rlbwt_row_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/rlbwt_structure_test: ./tests/unit/rlbwt/rlbwt_structure_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/permutation_lf_fl_test: ./tests/unit/rlbwt/permutation_lf_fl_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/permutation_phi_phi_inv_test: ./tests/unit/rlbwt/permutation_phi_phi_inv_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/rlbwt_interval_encoding_test: ./tests/unit/rlbwt/rlbwt_interval_encoding_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Integration-style tests that exercise larger rlbwt/runperm flows
$(INTEGRATION_BUILD_DIR)/rlbwt_test: ./tests/integration/rlbwt_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRATION_BUILD_DIR)/move_structure_test: ./tests/integration/move_structure_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRATION_BUILD_DIR)/move_test: ./tests/integration/move_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRATION_BUILD_DIR)/permutation_test: ./tests/integration/permutation_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Benchmarks (not run by default in `make test`)
$(BENCH_BUILD_DIR)/move_bench: ./tests/benchmark/move_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BENCH_BUILD_DIR)/permutation_bench: ./tests/benchmark/permutation_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BENCH_BUILD_DIR)/rlbwt_bench: ./tests/benchmark/rlbwt_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Root tests/benchmarks/example pick up -g -O0 via target-specific CXXFLAGS; ms/bi use their own Makefiles.
debug: CXXFLAGS += -g -O0
debug: all
	$(MAKE) -C src/ms debug
	$(MAKE) -C src/bi debug

clean:
	rm -rf $(BUILD_DIR)
	rm -f example
	$(MAKE) -C src/ms clean
	$(MAKE) -C src/bi clean
