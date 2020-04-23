CC=gcc
CFLAGS=-std=c11 -fno-strict-aliasing -msse2 -mfpmath=sse -D_GNU_SOURCE -g -Wall -I../../include
OPTFLAGS=-O3 -ffast-math -fno-finite-math-only -march=native
# OPTFLAGS=-O0 -g -ggdb3

all: obj2geo geo-subdiv geo-disp geo-make geo-createuvs obj2particles

obj2geo: obj2geo.c Makefile
	$(CC) -fopenmp $(CFLAGS) $(OPTFLAGS) ../../src/prims.c obj2geo.c -o obj2geo -lm

obj2particles: obj2particles.c Makefile
	$(CC) -fopenmp $(CFLAGS) $(OPTFLAGS) ../../src/prims.c obj2particles.c -o obj2particles -lm

geo-disp: disp.c texture.h vdata.h Makefile
	$(CC) -fopenmp $(CFLAGS) $(OPTFLAGS) ../../src/prims.c disp.c -o geo-disp -lm

geo-subdiv: subdiv.c subdiv.h vdata.h Makefile
	$(CC) -fopenmp $(CFLAGS) $(OPTFLAGS) ../../src/prims.c subdiv.c -o geo-subdiv -lm

geo-make: makegeo.c Makefile
	$(CC) -fopenmp $(CFLAGS) $(OPTFLAGS) makegeo.c -lm -o geo-make

geo-createuvs: createuvs.c Makefile
	$(CC) $(CFLAGS) $(OPTFLAGS) createuvs.c -lm -o geo-createuvs

clean:
	rm -f geo-disp geo-subdiv geo-make geo-createuvs obj2geo
