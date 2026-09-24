# MS application: build from project root with make -C src/ms
CXX = g++
CXXFLAGS = -std=c++17 -I../../include -I. -O3

MS_SRCS = ms.cpp tsv.cpp ms_io.cpp ms_test.cpp inspect.cpp tms_test.cpp

MS_HDRS = $(wildcard *.hpp) $(shell find ../../include -name '*.hpp')

ms: $(MS_SRCS) $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -o ms $(MS_SRCS)

# Debug build: -O0 -g for breakpoints (e.g. lldb ./ms, then run test ../../data)
debug: CXXFLAGS += -O0 -g
debug: ms

# Microbenchmark of LF and psi step chains on a tms index
tms_chain_bench: tms_chain_bench.cpp $(MS_HDRS)
	$(CXX) $(CXXFLAGS) -o tms_chain_bench tms_chain_bench.cpp

.PHONY: test all clean debug
test: ms
	./ms test ../../data

all: ms

clean:
	rm -f ms tms_chain_bench
