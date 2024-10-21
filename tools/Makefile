CC=clang
CFLAGS+=-Wall -std=c11 -msse2 -mfpmath=sse -fno-strict-aliasing -D_GNU_SOURCE -g -I. -I.. -I../build -I../include 
OPTFLAGS=-O3 -ffast-math -fno-finite-math-only -march=native

all: corona-netrender turntablecam battle-test
	make -C img/
	make -C geo/

corona-netrender: networkplayer.c ../src/display.d/xorg.c Makefile ../include/display.h
	$(CC) $(CFLAGS) $(OPTFLAGS) networkplayer.c ../src/display.d/xorg.c -lX11 -lrt -ljpeg -lm -lc -o corona-netrender

turntablecam: turntablecam.c Makefile ../include/camera.h ../include/quaternion.h
	$(CC) $(CFLAGS) $(OPTFLAGS) turntablecam.c -o turntablecam -lm

lerpcam: lerpcam.c Makefile ../include/camera.h ../include/quaternion.h
	$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@ -lm

# FIXME: only works with random numbers, not halton, only works with default PATHSPACE_MAX_VERTS etc..
battle-test: battle-test.c Makefile ../include/pathspace.h
	$(CC) $(CFLAGS) -I../ext/pthread-pool $(OPTFLAGS) battle-test.c ../src/pathspace.c ../src/lights.d/list.c ../src/camera.d/thinlens.c ../src/prims.c ../src/display.d/null.c ../src/render.d/gi.c ../src/sampler.d/pt.c ../src/pointsampler.d/rand.c ../src/points.d/sfmt.c ../src/shader.c ../src/view.c ../src/corona_common.c  -rdynamic -lm -lc -ldl -o battle-test -L../ext/pthread-pool/ -lpthreadpool

denoise: denoise.c denoise_wavelets.h Makefile
	$(CC) -std=c99 -g -I../include -lm -lc -Wall denoise.c -o denoise

atrous: atrous.c
	$(CC) $(OPTFLAGS) $(CFLAGS) atrous.c -o atrous -lm

fracf: fracf.c Makefile
	gcc -I.. -I../build -fno-strict-aliasing -fopenmp -march=native -mfpmath=sse -D_GNU_SOURCE -std=c11 -Wall -g -O3 -ffast-math fracf.c -lm -o fracf

glare: glare.c Makefile
	 gcc -fopenmp -fno-strict-aliasing -O3 -ffast-math -g -Wall -std=c99 -D_GNU_SOURCE -I ../build/ glare.c -lm -o glare

clean:
	rm -f corona-netrender battle-test turntablecam
