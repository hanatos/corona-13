CXX=clang++
CC=clang
INC=-I../../build -I../../include -I../.. -I.
CFLAGS=-std=c11 -msse2 -mfpmath=sse -D_GNU_SOURCE -g -Wall -fno-strict-aliasing $(INC)
CXXFLAGS=-msse2 -mfpmath=sse -D_GNU_SOURCE -g -Wall -fno-strict-aliasing
OPTFLAGS=-O3 -ffast-math -fno-finite-math-only -march=native
EXRINC=-I/usr/include/OpenEXR
EXRLD=-lIex -lIlmImf -lImath -lHalf -lz

all: pfm2exr pfm2hdr bump2normalmap exr2dng exr2pfm pfmdiff bakesky welch pfm2fb ppm2fb rgb2spec_opt

rgb2spec_opt: rgb2spec_opt.cpp ../../include/rgb2spec.h details/cie1931.h details/lu.h
	$(CXX) -std=c++11 -fopenmp $< -o $@ -lomp 

welch: welch.c Makefile
	$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@ -lm

pfm2exr: pfm2exr.cpp Makefile
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(EXRINC) pfm2exr.cpp -o pfm2exr $(EXRLD)

# disabled now, needs to be discarded or rewritten (see exr2dng)
# pfm2dng: pfm2dng.c ../../include/screenshot_dng.h Makefile
# $(CC) $(CFLAGS) $(OPTFLAGS) pfm2dng.c -o pfm2dng -lm

pfm2hdr: pfm2hdr.c rgbe.c rgbe.h Makefile
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INC) rgbe.c pfm2hdr.c -o pfm2hdr -lm

pfmdiff: pfmdiff.c Makefile
	$(CC) $(CFLAGS) $(OPTFLAGS) pfmdiff.c -o pfmdiff -lm

bump2normalmap: bump2normalmap.c Makefile
	$(CC) $(CFLAGS) $(OPTFLAGS) bump2normalmap.c -o bump2normalmap -lm

exr2dng: exr2dng.cpp Makefile
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(EXRINC) exr2dng.cpp -o exr2dng $(EXRLD)

exr2pfm: exr2pfm.cpp Makefile
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(EXRINC) exr2pfm.cpp -o exr2pfm $(EXRLD)

# needs rewrite/discard
# exr2hdr: exr2hdr.cpp rgbe.c rgbe.h Makefile
# $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(EXRINC) exr2hdr.cpp -o exr2hdr $(EXRLD)

# needs rewrite/discard
# spec2rgb: spec2rgb.c rgbe.c rgbe.h Makefile
# 	$(CC) $(CFLAGS) $(OPTFLAGS) $(INC) rgbe.c spec2rgb.c -o spec2rgb

dmap2pgm: dmap2pgm.c Makefile
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INC) dmap2pgm.c -o dmap2pgm

bakesky: bakesky.c Makefile
	 $(CC) -fopenmp $(CFLAGS) $(OPTFLAGS) $(INC) bakesky.c ../../src/screenshot.c -o bakesky -lm

dbor: dbor.c Makefile
	 $(CC) $(CFLAGS) $(OPTFLAGS) $(INC) $< -o $@ -lm

ppm2fb: ppm2fb.c Makefile
	 $(CC) $(CFLAGS) $(OPTFLAGS) $(INC) $< -o $@ -lm

dbord: dbord.c Makefile
	 $(CC) $(CFLAGS) -fopenmp $(OPTFLAGS) $(INC) $< -o $@ -lm -lomp

clean:
	rm -f pfm2exr pfm2dng pfm2hdr pfmdiff bump2normalmap exr2dng exr2pfm exr2hdr spec2rgb dmap2pgm bakesky dbor
