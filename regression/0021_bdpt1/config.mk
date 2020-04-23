
# configure your favorite rt here.

# kdtree, bih, qbvh, grid
MOD_accel=qbvhmp
# thinlens, lens
MOD_camera=thinlens
# LDFLAGS+=../optics/gui-backend.o ../optics/liblensflare.a -lstdc++
# CFLAGS+=-I. -I../optics
# in case you do realistic lenses, choose the one here.
# see subdirectories of camera/ for a list.
CAMERA_LENS=wideangle-1971
# null, network, unix, gl
MOD_display=null
# LDFLAGS+=-lX11 # for unix display
# LDFLAGS+=-lGL -lSDL # for gl
# LDFLAGS+=-lswscale -lavcodec # for network display
# sfmt, sobol
MOD_points=sfmt
# LDFLAGS+=-lz # for sobol points
# rdtsc, vis, gi, tiles
MOD_render=gi
# pixel filter: box, bilin, spline
BUILD_MOD_filter=blackmanharris
# pt, ptdl, lt, ptlt, bdpt, ppm
MOD_sampler=bdpt1
# screenshot format: pfm, dng
MOD_screenshot=pfm
# colorspace of camera: xyz, rgb
COL_camera=xyz
COL_input=ergb
COL_output=srgb
MOD_spectrum=grid
# point sampling: rand, kmlt
MOD_pointsampler=rand
# space separated list of mutation strategies
MUTATIONS=largestep halfvec lens

include arch
