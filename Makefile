CC=gcc
CFLAGS+=-fPIC -fno-strict-aliasing -std=c11 -Wall -pipe -I/usr/include -Ibuild/ -I. -Iinclude/ -Icamera/ -Iext/pthread-pool/ -D_GNU_SOURCE -g
# openvdb will be appended for all but -I
OPENVDB_PATH=ext/openvdb_dev
#CFLAGS+=-Werror -Wno-error=maybe-uninitialized
LDFLAGS=-lm -lc -ldl -rdynamic -Lext/pthread-pool/ -lpthreadpool -pthread -Wl,-rpath,'$$ORIGIN/shaders'
SHARED=-shared
# dr dobb's idea about makefile debugging:
OLD_SHELL := $(SHELL)
# SHELL = $(warning [$@ ($^) ($?)])$(OLD_SHELL)
SHELL = $(warning [$@ ($?)])$(OLD_SHELL)

.PHONY: all clean debug modules modules_clean

include config.mk
CFLAGS+=-Icamera/$(CAMERA_LENS)/
CFLAGS+=-Isrc/pointsampler.d/ # for vmlt.h

MOD_HEADERS=\
    build/filter.h \
    build/vmlt_registry.h\
    build/colourspaces.h

MOD_SOURCES=\
    $(EXTRA_SOURCES)\
    src/accel.d/$(MOD_accel).c\
    src/camera.d/$(MOD_camera).c\
    src/display.d/$(MOD_display).c\
    src/points.d/$(MOD_points).c\
    src/pointsampler.d/$(MOD_pointsampler).c\
    src/render.d/$(MOD_render).c\
    src/sampler.d/$(MOD_sampler).c\
    src/lights.d/$(MOD_lights).c

CAMERA_HEADERS=\
    camera/$(CAMERA_LENS)/init.h\
    camera/$(CAMERA_LENS)/pt_evaluate.h\
    camera/$(CAMERA_LENS)/pt_evaluate_aperture.h\
    camera/$(CAMERA_LENS)/pt_sample_aperture.h\
    camera/$(CAMERA_LENS)/lt_sample_aperture.h\
    camera/$(CAMERA_LENS)/pt_evaluate_jacobian.h\
    camera/$(CAMERA_LENS)/pt_evaluate_aperture_jacobian.h

STATIC_HEADERS=\
    camera/lens.h\
    include/corona_common.h\
    include/camera.h\
    include/render.h\
    include/view.h\
    include/points.h\
    include/sampler.h\
    include/spectrum.h\
    include/pointsampler.h\
    include/lbvh.h\
    include/dbor.h\
    include/fakegaussian.h\
    include/distancemap.h\
    include/quaternion.h\
    include/rgb2spec.h\
    include/shader.h\
    include/prims.h \
    include/geo.h \
    include/geo/line.h\
    include/geo/sphere.h\
    include/geo/triangle.h\
    include/pathspace/manifold.h\
    include/pathspace/measurement.h\
    include/pathspace/raydifferentials.h\
    include/pathspace/nee.h\
    include/pathspace/mnee.h\
    include/pathspace/tech.h\
    include/screenshot.h\
    include/lights.h

STATIC_SOURCES=\
    src/corona_common.c\
    src/shader.c\
    src/view.c\
    src/pathspace.c\
    src/pathspace/tech.c\
    src/prims.c\
    src/screenshot.c\
    src/main.c

# veach mlt subsystem required for vmlt and erpt
ifeq ($(MOD_pointsampler),vmlt)
	STATIC_SOURCES+=src/pathspace/vmlt.c
endif
ifeq ($(MOD_render),erpt)
	STATIC_SOURCES+=src/pathspace/vmlt.c
endif

debug:CFLAGS+=-gdwarf-2 -ggdb3 -msse2 -mfpmath=sse -O0
debug:LDFLAGS+=-lm -lc -ldl -lrt -lX11 -rdynamic -Lext/pthread-pool/ -lpthreadpool -pthread -Wl,-rpath,'$$ORIGIN/shaders'
debug:corona modules

# use address sanitizer feature of clang. run with:
# ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./corona
# or whatever the path in your system if it's not in $PATH.
sanitize:CFLAGS+=-O0 -ggdb3 -msse2 -mfpmath=sse -fsanitize=address -fno-omit-frame-pointer 
# sanitize:CFLAGS+=-O1 -fsanitize=address -fno-omit-frame-pointer
sanitize:LDFLAGS+=-fsanitize=address
sanitize:corona modules

config.mk arch:
	@echo -e " \033[1m\033[33m*\033[30m\033[0m please configure the build for your processor and rendering needs:"
	@echo -e " \033[1m\033[33m*\033[30m\033[0m cp arch.example arch; cp config.mk.example config.mk" && false

build/%.h : include/%*.h Makefile config.mk
	@cp -f include/$*_$(value BUILD_MOD_$*).h $@
build/%.c : src/%*.c Makefile config.mk
	@cp -f src/$*_$(value BUILD_MOD_$*).c $@

# configure colour spaces
build/colourspaces.h: Makefile config.mk
	@echo "#pragma once" > build/colourspaces.h
	@echo "// generated during build, do not edit." >> build/colourspaces.h
	@for file in include/colour/*.h; do \
		echo \#include \"$$file\" >> build/colourspaces.h;\
	done
	@echo "#define colour_input_to_xyz colour_${COL_input}_to_xyz" >> build/colourspaces.h
	@echo "#define colour_xyz_to_input colour_xyz_to_${COL_input}" >> build/colourspaces.h
	@echo "#define colour_input_print_info colour_${COL_input}_print_info" >> build/colourspaces.h
	@echo "#define colour_output_to_xyz colour_${COL_output}_to_xyz" >> build/colourspaces.h
	@echo "#define colour_xyz_to_output colour_xyz_to_${COL_output}" >> build/colourspaces.h
	@echo "#define colour_output_print_info colour_${COL_output}_print_info" >> build/colourspaces.h
	@echo "#define colour_camera_to_xyz colour_${COL_camera}_to_xyz" >> build/colourspaces.h
	@echo "#define colour_xyz_to_camera colour_xyz_to_${COL_camera}" >> build/colourspaces.h
	@echo "#define colour_camera_print_info colour_${COL_camera}_print_info" >> build/colourspaces.h

build:
	mkdir -p build

# version header
build/version.h: Makefile .git/FETCH_HEAD
	@echo "#ifndef CORONA_VERSION_H" > build/version.h
	@echo "#define VERSION \"$(shell git describe --tags)\"" >> build/version.h
	@echo "#endif" >> build/version.h

# veach metropolis configuration
build/vmlt_registry.h: Makefile config.mk
	@printf "%s\n%s\n%s\n%s\n%s\n%s\n%s\n" $(foreach MUT,$(MUTATIONS),"#include \"vmlt_$(MUT).h\"") > build/vmlt_registry.h
	@printf "static inline void vmlt_register_all(vmlt_t *t)\n{\n" >> build/vmlt_registry.h
	@printf "  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n" $(foreach MUT,$(MUTATIONS),"vmlt_register(t, &$(MUT)_init, &$(MUT)_cleanup, &$(MUT)_suitability, &$(MUT)_mutate, &$(MUT)_print_info);") >> build/vmlt_registry.h
	@echo "}" >> build/vmlt_registry.h

build/filter.h: Makefile config.mk
	@printf "#pragma once\n" > build/filter.h
	@printf '#include "filter/box.h"\n' >> build/filter.h
	@printf '#include "filter/blackmanharris.h"\n' >> build/filter.h
	@printf "#define filter_splat filter_%s_splat\n" $(MOD_filter) >> build/filter.h
	@printf "#define filter_splat4 filter_%s_splat4\n" $(MOD_filter) >> build/filter.h
	@printf "#define filter_print_info filter_%s_print_info\n" $(MOD_filter) >> build/filter.h

ext/pthread-pool/libpthreadpool.a:
	make -C ext/pthread-pool/

test:
	make -C regression

# create tristimulus to spectrum upsampling luts:
COLOURLUTS=data/ergb2spec.coeff data/xyz2spec.coeff
tools/img/rgb2spec_opt:
	make -C tools/img rgb2spec_opt

data/ergb2spec.coeff: tools/img/rgb2spec_opt
	mkdir -p data
	tools/img/rgb2spec_opt 64 $@ eRGB

data/xyz2spec.coeff: tools/img/rgb2spec_opt
	mkdir -p data
	tools/img/rgb2spec_opt 64 $@ XYZ

corona: build $(MOD_HEADERS) $(MOD_SOURCES) $(COLOURLUTS) $(CAMERA_HEADERS) $(STATIC_HEADERS) $(STATIC_SOURCES) Makefile arch build/version.h ext/pthread-pool/libpthreadpool.a
	$(CC) $(CFLAGS) $(STATIC_SOURCES) $(MOD_SOURCES) -o corona $(LDFLAGS)

clean: modules_clean
	rm -f corona
	rm -f build/*
	make -C ext/pthread-pool clean

SHADERS=color mult interior texture medium_rgb medium_poe dielectric diffdiel colorcheckersg metal medium_hete vdata medium_aggregate bump
ifneq (,$(findstring MF_COUNT,$(CFLAGS)))
  # hero wavelength, need to disable a few non-ported bsdf :(
else
  # single wavelength
  SHADERS+=mmetal mdiffuse mdielectric hair
endif
SHADER_DEPEND=include/corona_common.h include/spectrum.h build/colourspaces.h include/shader.h Makefile
shaders/libdielectric.so:src/shaders/ggx.h
shaders/libmdielectric2.so:src/shaders/MicrosurfaceScattering.cpp src/shaders/MicrosurfaceScattering.h
shaders/libmmetal2.so:src/shaders/MicrosurfaceScattering.cpp src/shaders/MicrosurfaceScattering.h
shaders/libmdiffuse2.so:src/shaders/MicrosurfaceScattering.cpp src/shaders/MicrosurfaceScattering.h
shaders/libdiffdiel.so:src/shaders/ggx.h
shaders/libmdielectric.so:src/shaders/microfacet.h
shaders/libmetal.so:src/shaders/fresnel.h src/shaders/ggx.h
shaders/libmmetal.so:src/shaders/fresnel.h src/shaders/microfacet.h
shaders/libmdiffuse.so:src/shaders/microfacet.h

VOL_HEADERS=include/vol/types.h\
            include/vol/vol.h\
            include/vol/trace.h\
            include/vol/trace_impl.inc\
            include/vol/interpolation.h\
            include/vol/payload.h\
            include/vol/payload_compress.h\
            include/vol/lighthierarchy.h\
            include/vol/shaders.h\
            include/vol/trace_octree.h
shaders/libmedium_hete.so:${VOL_HEADERS}

modules_clean:
	rm -f shaders/*.so

modules: build $(MOD_HEADERS) $(MOD_SOURCES)
modules: CFLAGS+=-fPIC
modules: $(patsubst %,shaders/lib%.so,$(SHADERS)) \
	shaders/libenvmap.so \
  shaders/libconst.so
#shaders/libspectral.so  # this is cool and we want it back!

shaders/lib%.so: src/shaders/%.c $(SHADER_DEPEND) # shaders
	@mkdir -p shaders
	$(CC) $(CFLAGS) $(SHARED) $< -o $@

shaders/libenvmap.so: src/shaders/sky_envmap.c $(SHADER_DEPEND)
	$(CC) $(CFLAGS) $(SHARED) $< -I src/shaders/ -o $@
shaders/libspectral.so: src/shaders/sky_spectral.c $(SHADER_DEPEND)
	$(CC) $(CFLAGS) $(SHARED) $< -I src/shaders/ -o $@
shaders/libconst.so: src/shaders/sky_const.c $(SHADER_DEPEND)
	$(CC) $(CFLAGS) $(SHARED) $< -o $@


