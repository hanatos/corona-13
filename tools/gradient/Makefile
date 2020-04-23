
.PHONY:all clean

all:poissonsolve subgrad cgsolve

CPP=clang++
CC=clang-6.0
CFLAGS=-std=c99 -Wall -march=native -fopenmp -O3
CPPFLAGS=-std=c++11 -Wall -march=native -I/usr/include/eigen3/ -g
LDFLAGS=-lm -lfftw3f -lfftw3f_threads

poissonsolve: main.c fft_solver.h Makefile
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

cgsolve: main.c cgsolve.h Makefile
	$(CC) $(CFLAGS) -fsave-optimization-record $< -o $@ -lm

hubersolve: hubersolve.cc hubersolve.hh Makefile
	$(CPP) $(CPPFLAGS) $< -o $@

subgrad: subgrad.c Makefile
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f poissonsolve subgrad
