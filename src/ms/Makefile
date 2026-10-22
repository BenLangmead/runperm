# MS application: build from project root with make -C src/ms
CXX = g++
CXXFLAGS = -std=c++17 -I../../include -I. -O3 -pthread

# The commit recorded in the ftabs this binary builds; outside a git checkout
# (as from git archive), pass MS_COMMIT=... on the command line.
MS_COMMIT ?= $(shell c=$$(git rev-parse --short=12 HEAD 2>/dev/null) && { git diff --quiet HEAD -- . 2>/dev/null || c=$$c-dirty; echo $$c; } || echo unknown)
CXXFLAGS += -DMS_COMMIT='"$(MS_COMMIT)"'

MS_SRCS = ms.cpp ms_io.cpp index_format.cpp ms_test.cpp inspect.cpp tms_test.cpp fmd_test.cpp

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

# Benchmark build that also records, for ms batch and tms-batch with
# --interleave 0, where in each read the repositions, the rows they read and
# the query time fall (see step_stats.hpp), as ms-steps.
ms-steps: $(MS_SRCS) $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DMS_STATS -DTMS_STATS -DMS_STEP_STATS -o ms-steps $(MS_SRCS) -lz

# Debug build: -O0 -g for breakpoints (e.g. lldb ./ms, then run test ../../data)
debug: CXXFLAGS += -O0 -g
debug: ms

# Microbenchmark of LF and psi step chains on a tms index
tms_chain_bench: tms_chain_bench.cpp index_format.cpp $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -o tms_chain_bench tms_chain_bench.cpp index_format.cpp

# Benchmark of smem-all's counting methods on one load of a tms index (see
# tms_count_bench.cpp); counts matching statistics repositions as ms-steps does
tms_count_bench: tms_count_bench.cpp index_format.cpp $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DTMS_STATS -o tms_count_bench tms_count_bench.cpp index_format.cpp -lz

.PHONY: test test-threads test-reader test-rewrap test-ftab all clean debug bench
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

# Stored ftabs: ftab-build, ftab-verify and byte identity with --ftab-file
test-ftab: ms
	./test_ftab.sh ./ms ../../data

# index-rewrap against new builds: OLD is an ms binary that writes the
# earlier formats (MSIX v4, TMSX v8), as built from commit f82f172
test-rewrap: ms
	@test -n "$(OLD)" || { echo "usage: make test-rewrap OLD=path/to/old/ms"; exit 1; }
	./test_rewrap.sh $(OLD) ./ms ../../data

# tms-batch --pos-format rope against ropebwt3 mem on a small FMD: TERALCP
# is a TeraLCP that writes -ophi, and ROPEBWT3 the ropebwt3 to compare with
test-rope: ms tools/fmd2rlbwt/fmd2rlbwt
	@test -n "$(TERALCP)" || { echo "usage: make test-rope TERALCP=path/to/TeraLCP [ROPEBWT3=ropebwt3]"; exit 1; }
	./test_rope.sh ./ms tools/fmd2rlbwt/fmd2rlbwt $(TERALCP) $(or $(ROPEBWT3),ropebwt3)

tools/fmd2rlbwt/fmd2rlbwt:
	$(MAKE) -C tools/fmd2rlbwt

all: ms

clean:
	rm -f ms ms-bench ms-stats ms-steps tms_chain_bench read_chunks_test
