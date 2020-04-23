CXX=g++
CC=gcc
CFLAGS+=-std=c11 -msse2 -mfpmath=sse -D_GNU_SOURCE -g -Wall -fno-strict-aliasing
CXXFLAGS=-std=c++11 -msse2 -mfpmath=sse -D_GNU_SOURCE -g -Wall -fno-strict-aliasing
LDFLAGS=-lX11 -lrt -ljpeg
OPTFLAGS=-O3 -fopenmp -ffast-math -fno-finite-math-only -march=native # -DNDEBUG
# OPTFLAGS=-gdwarf-2 -g3
# CFLAGS+=-fsave-optimization-record
# CFLAGS+=-pg
# LDFLAGS+=-pg

all: voltest vollh slice

VOLHEADERS=../../include/vol/types.h ../../include/vol/trace.h ../../include/vol/trace_impl.inc ../../include/vol/vol.h ../../include/vol/interpolation.h ../../include/vol/payload.h ../../include/vol/payload_compress.h ../../include/vol/lighthierarchy.h
SPEC=../../src/spectrum.d/grid.c

voltest: voltest.c ${VOLHEADERS} Makefile ../../include/points.h cloud.h
	$(CC) $(OPTFLAGS) $(CFLAGS) -I. -I../.. -I ../../include/ -I../../build/ $< $(SPEC) ../../src/points.d/sfmt.c -lm -o $@ 

slice: slice.c ${VOLHEADERS} Makefile
	$(CC) $(OPTFLAGS) $(CFLAGS) -I. -I../.. -I ../../include/ -I../../build/ $< $(SPEC) -lm -o $@

vollh: vollh.c ${VOLHEADERS} Makefile ../../include/points.h
	$(CC) $(OPTFLAGS) $(CFLAGS) -I. -I../.. -I ../../include/ -I../../build/ $< $(SPEC) ../../src/points.d/sfmt.c -lm -o $@

ptc2vol: ptc2vol.c ${VOLHEADERS} Makefile
	clang -fno-strict-aliasing -g $(OPTFLAGS) -Wall -std=c11 -D_GNU_SOURCE -I. -I ../../include/ -I../../build/ ptc2vol.c -lm -o ptc2vol

OPENVDB_LIB=/usr/lib/x86_64-linux-gnu
OPENVDB_INC=/usr/include/
ovdb2vol: ovdb2vol.cc ${VOLHEADERS} Makefile
	 $(CXX) -fopenmp ${CXXFLAGS} ${OPTFLAGS} -I../../build -I../../include/ -I../../ -I$(OPENVDB_INC) ovdb2vol.cc $(SPEC) -lm -L$(OPENVDB_LIB) -lopenvdb -ltbb -lboost_system -lboost_iostreams -lHalf -Wl,-rpath,'$$ORIGIN/$(OPENVDB_LIB)' -o ovdb2vol

clean:
	rm -f voltest ovdb2vol ptc2vol slice vollh
