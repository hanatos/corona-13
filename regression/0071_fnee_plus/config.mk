# qbvh, qbvhm (motion blur), qbvhmp (motion + parallel build)
MOD_accel=qbvhmp
# thinlens, polynomial
MOD_camera=thinlens
# in case you do realistic lenses, choose the one here.
# see subdirectories of camera/ for a list.
CAMERA_LENS=wideangle-1971
# null, mjpeg, unix, gl
MOD_display=null
# LDFLAGS+=-lX11 # for xorg display
# LDFLAGS+=-lGL -lSDL # for gl
# LDFLAGS+=-ljpeg # mjpeg network display
# sfmt, xorshift128p
MOD_points=sfmt
# LDFLAGS+=-lz # for sobol points
# vis, gi, tiles
MOD_render=gi
# pixel filter: box, bilin, spline
BUILD_MOD_filter=blackmanharris
# pt, ptdl, lt, ptlt, bdpt, ppm
MOD_sampler=ptdl
# rgb to spectrum upsampling method: smits, grid
MOD_spectrum=grid
# point sampling: rand, kmlt, vmlt, halton
MOD_pointsampler=rand
# space separated list of mutation strategies for vmlt
MUTATIONS=largestep multichain
# colour configuration
# input can be anything: ergb (for smits), xyz, aces, srgb, rec709
COL_input=ergb
# output can have gamma: adobergb, srgb, or custom profile
COL_output=adobergb
# framebuffer needs to be linear: xyz, rec709
COL_camera=xyz

include arch
# if overwriting this, be sure to `make clean':
# CFLAGS+=-DNDEBUG
# CFLAGS+=-DPATHSPACE_MAX_VERTS=3
# CXXFLAGS+=-DPATHSPACE_MAX_VERTS=2

# stupid regularisation based on path length:
# CFLAGS+=-DSHADER_ROUGHENING
CFLAGS+=-DFNEE -DSEGMENT_EMISSION
