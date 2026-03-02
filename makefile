VERSION = 0.1.3

CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O3

HEADERS = $(shell find include -name "*.hpp")

BUILD_DIR = build
UNIT_BUILD_DIR = $(BUILD_DIR)/unit
INTEGRATION_BUILD_DIR = $(BUILD_DIR)/integration
BENCH_BUILD_DIR = $(BUILD_DIR)/bench

# High-level targets
# - all: build all unit/integration tests + benchmarks + examples
# - test: build unit/integration tests and run them

UNIT_TESTS = $(UNIT_BUILD_DIR)/packed_vector_test
INTEGRATION_TESTS = $(INTEGRATION_BUILD_DIR)/rlbwt_test
BENCH_TESTS = $(BENCH_BUILD_DIR)/move_bench \
              $(BENCH_BUILD_DIR)/runperm_bench \
              $(BENCH_BUILD_DIR)/rlbwt_bench

# all: build invert move_build
all: $(UNIT_TESTS) $(INTEGRATION_TESTS) $(BENCH_TESTS) example_test

test: $(UNIT_TESTS) $(INTEGRATION_TESTS)
	$(UNIT_BUILD_DIR)/packed_vector_test
	$(INTEGRATION_BUILD_DIR)/rlbwt_test

example_test: examples/example_test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Unit tests (header-only data structures)
$(UNIT_BUILD_DIR)/packed_vector_test: ./tests/unit/ds/packed_vector_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Integration-style tests that exercise larger rlbwt/runperm flows
$(INTEGRATION_BUILD_DIR)/rlbwt_test: ./tests/integration/rlbwt_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Benchmarks (not run by default in `make test`)
$(BENCH_BUILD_DIR)/move_bench: ./tests/benchmark/move_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BENCH_BUILD_DIR)/runperm_bench: ./tests/benchmark/runperm_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BENCH_BUILD_DIR)/rlbwt_bench: ./tests/benchmark/rlbwt_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY: debug
debug: CXXFLAGS += -g -O0
debug: all

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f invert move_build example_test