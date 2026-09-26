# MS application: build from project root with make -C src/ms
CXX = g++
CXXFLAGS = -std=c++17 -I../../include -I. -O3 -pthread

MS_SRCS = ms.cpp tsv.cpp ms_io.cpp ms_test.cpp inspect.cpp tms_test.cpp

MS_HDRS = $(wildcard *.hpp) $(shell find ../../include -name '*.hpp')

ms: $(MS_SRCS) $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -o ms $(MS_SRCS) -lz

# Benchmark build: assertions off, as ms-bench.  The test commands refuse to
# run in it, since the suites check their results with assert.
ms-bench: $(MS_SRCS) $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -o ms-bench $(MS_SRCS) -lz

bench: ms-bench

# Benchmark build that also counts LF fast-forward rows and reposition walk
# lengths in ms batch (see MsStats), as ms-stats.
ms-stats: $(MS_SRCS) $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DMS_STATS -o ms-stats $(MS_SRCS) -lz

# Debug build: -O0 -g for breakpoints (e.g. lldb ./ms, then run test ../../data)
debug: CXXFLAGS += -O0 -g
debug: ms

# Microbenchmark of LF and psi step chains on a tms index
tms_chain_bench: tms_chain_bench.cpp $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -o tms_chain_bench tms_chain_bench.cpp

.PHONY: test test-threads test-reader all clean debug bench
test: ms
	./ms test ../../data

# Read chunking at every segment boundary, under the address and undefined
# behavior sanitizers
read_chunks_test: read_chunks_test.cpp read_chunks.hpp
	$(CXX) -std=c++17 -I. -O1 -g -pthread -fsanitize=address,undefined -fno-sanitize-recover=all -o read_chunks_test read_chunks_test.cpp

test-reader: read_chunks_test
	./read_chunks_test

# Output identity across interleave, thread count and block size
test-threads: ms
	./test_threads.sh ./ms ../../data

all: ms

clean:
	rm -f ms ms-bench ms-stats tms_chain_bench read_chunks_test
