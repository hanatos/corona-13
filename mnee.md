# instructions for testing manifold next event estimation (MNEE)

## compile
compile with `config.mk` maybe similar to:

```
MOD_pointsampler?=halton
MOD_sampler?=ptmnee
CFLAGS+=-DPATHSPACE_MAX_VERTS=10
```

you can use the xorg output module if you want an interactive window, in which
case you need to press `p` to print an image to an output buffer.

if you want to run command line, set

```
MOD_display?=null
```

instead of `xorg`. you can leave away `-lX11` then of course.

## run
because MNEE is not a complete technique, there will be fireflies
left in the render, caused by pt because ptmnee does not cover these paths.
to remove these, you'll want to run

```
./corona scenes/tumbler/tumbler.nra2 --dbor 11
```

or similar, to enable a couple of dbor layers.

to arrive at the firefly removed image, use

```
./tools/img/dbor scenes/tumbler/tumblerrender 0.02 -100
```

the output pfm is in XYZ colour space, look at it with whatever software
you like. the minimal semi-useful viewer is probably [`eu`](https://github.com/hanatos/eu).

you can repeat the above steps with `MOD_sampler?=ptdl` to compare against path tracing
with vanilla next event estimation only.
