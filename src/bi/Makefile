# BI application: build from project root with make -C src/bi
CXX = g++
CXXFLAGS = -std=c++17 -I../../include -I. -I../ms -O3

BI_SRCS = bi.cpp bi_test.cpp bi_io.cpp ../ms/tsv.cpp

bi: $(BI_SRCS)
	$(CXX) $(CXXFLAGS) -o bi $(BI_SRCS)

# Debug build: -O0 -g for breakpoints (e.g. lldb ./bi, then run test)
debug: CXXFLAGS += -O0 -g
debug: bi

.PHONY: test all clean debug
test: bi
	./bi test

all: bi

clean:
	rm -f bi
