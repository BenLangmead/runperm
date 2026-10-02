# fmd2rlbwt: dump a ropebwt3 FMD's maximal runs (see FMD_RUNS.md).
# rld0.c and rld0.h are vendored from lh3/ropebwt3 (MIT, LICENSE.ropebwt3).
CC = gcc
CXX = g++
CFLAGS = -O2 -Wall
CXXFLAGS = -std=c++17 -O2 -Wall

fmd2rlbwt: fmd2rlbwt.cpp rld0.o rld0.h
	$(CXX) $(CXXFLAGS) -o $@ fmd2rlbwt.cpp rld0.o

rld0.o: rld0.c rld0.h
	$(CC) $(CFLAGS) -c rld0.c

.PHONY: clean
clean:
	rm -f fmd2rlbwt rld0.o
