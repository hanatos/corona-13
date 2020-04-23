#!/bin/bash

# iccdump outputs the transpose of the rgb -> xyz matrix, so read the transpose of that:
read -r w0 w1 w2 a00 a10 a20 a01 a11 a21 a02 a12 a22 <<< $(iccdump -v3 -t wtpt -t rXYZ -t gXYZ -t bXYZ $1 | grep 0: | awk '{print $2 $3 $4}' | tr ',' ' ')
echo "wp $w0 $w1 $w2"
echo "rgb to xyz"
echo $a00 $a01 $a02
echo $a10 $a11 $a12
echo $a20 $a21 $a22
# get same chromaticities relative to illuminant E
# w0=1.0
# w1=1.0
# w2=1.0
# DEBUG: see adaptation matrix:
# a00=1.0
# a01=0.0
# a02=0.0
# a10=0.0
# a11=1.0
# a12=0.0
# a20=0.0
# a21=0.0
# a22=1.0
# bradford adaptation Ma
ma00=0.8951000;  ma01=0.2664000;  ma02=-0.1614000
ma10=-0.7502000; ma11=1.7135000;  ma12=0.0367000
ma20=0.0389000;  ma21=-0.0685000; ma22=1.0296000
# bradford adaptation Ma_inv
mai00=0.9869929;  mai01=-0.1470543; mai02=0.1599627
mai10=0.4323053;  mai11=0.5183603;  mai12=0.0492912
mai20=-0.0085287; mai21=0.0400428;  mai22=0.9684867

# create matrix m = mai s ma a:

# b = ma * a
b00=$(echo $ma00*$a00 + $ma01*$a10 + $ma02*$a20 | bc -l)
b01=$(echo $ma00*$a01 + $ma01*$a11 + $ma02*$a21 | bc -l)
b02=$(echo $ma00*$a02 + $ma01*$a12 + $ma02*$a22 | bc -l)
b10=$(echo $ma10*$a00 + $ma11*$a10 + $ma12*$a20 | bc -l)
b11=$(echo $ma10*$a01 + $ma11*$a11 + $ma12*$a21 | bc -l)
b12=$(echo $ma10*$a02 + $ma11*$a12 + $ma12*$a22 | bc -l)
b20=$(echo $ma20*$a00 + $ma21*$a10 + $ma22*$a20 | bc -l)
b21=$(echo $ma20*$a01 + $ma21*$a11 + $ma22*$a21 | bc -l)
b22=$(echo $ma20*$a02 + $ma21*$a12 + $ma22*$a22 | bc -l)

#  D50   0.96422   1.00000   0.82521
#  D55   0.95682   1.00000   0.92149
#  D65   0.95047   1.00000   1.08883
#  E     1.0       1.0       1.0
# icc profiles store their colorants relative to D50, correct for that:
x=0.96422; y=1.00000; z=0.82521
# x=1.0 y=1.0 z=1.0
d50l=$(echo $ma00*$x + $ma01*$y + $ma02*$z | bc -l)
d50m=$(echo $ma10*$x + $ma11*$y + $ma12*$z | bc -l)
d50s=$(echo $ma20*$x + $ma21*$y + $ma22*$z | bc -l)
wl=$(echo $ma00*$w0 + $ma01*$w1 + $ma02*$w2 | bc -l)
wm=$(echo $ma10*$w0 + $ma11*$w1 + $ma12*$w2 | bc -l)
ws=$(echo $ma20*$w0 + $ma21*$w1 + $ma22*$w2 | bc -l)
s0=$(echo $wl/$d50l | bc -l)
s1=$(echo $wm/$d50m | bc -l)
s2=$(echo $ws/$d50s | bc -l)
# c = s * b, where s is the scale matrix with the white point
c00=$(echo $s0*$b00 | bc -l)
c01=$(echo $s0*$b01 | bc -l)
c02=$(echo $s0*$b02 | bc -l)
c10=$(echo $s1*$b10 | bc -l)
c11=$(echo $s1*$b11 | bc -l)
c12=$(echo $s1*$b12 | bc -l)
c20=$(echo $s2*$b20 | bc -l)
c21=$(echo $s2*$b21 | bc -l)
c22=$(echo $s2*$b22 | bc -l)
# 
# m = mai * c
m00=$(echo $mai00*$c00 + $mai01*$c10 + $mai02*$c20 | bc -l)
m01=$(echo $mai00*$c01 + $mai01*$c11 + $mai02*$c21 | bc -l)
m02=$(echo $mai00*$c02 + $mai01*$c12 + $mai02*$c22 | bc -l)
m10=$(echo $mai10*$c00 + $mai11*$c10 + $mai12*$c20 | bc -l)
m11=$(echo $mai10*$c01 + $mai11*$c11 + $mai12*$c21 | bc -l)
m12=$(echo $mai10*$c02 + $mai11*$c12 + $mai12*$c22 | bc -l)
m20=$(echo $mai20*$c00 + $mai21*$c10 + $mai22*$c20 | bc -l)
m21=$(echo $mai20*$c01 + $mai21*$c11 + $mai22*$c21 | bc -l)
m22=$(echo $mai20*$c02 + $mai21*$c12 + $mai22*$c22 | bc -l)

#debug: with a=identity, this results in the same matrix as bruce lindbloom has as d50->d65 chromatic adaptation (bradford)
echo "bradford adapted rgb to xyz"
echo $m00 $m01 $m02
echo $m10 $m11 $m12
echo $m20 $m21 $m22

# now invert the matrix to get xyz -> rgb as needed in the code
invdet=$(echo 1/\($m00*\($m22*$m11 - $m21*$m12\) - $m10*\($m22*$m01 - $m21*$m02\) + $m20*\($m12*$m01 - $m11*$m02\) \)| bc -l);
i00=$(echo  $invdet*\($m22*$m11 - $m21*$m12\) | bc -l)
i01=$(echo -$invdet*\($m22*$m01 - $m21*$m02\) | bc -l)
i02=$(echo  $invdet*\($m12*$m01 - $m11*$m02\) | bc -l)
i10=$(echo -$invdet*\($m22*$m10 - $m20*$m12\) | bc -l)
i11=$(echo  $invdet*\($m22*$m00 - $m20*$m02\) | bc -l)
i12=$(echo -$invdet*\($m12*$m00 - $m10*$m02\) | bc -l)
i20=$(echo  $invdet*\($m21*$m10 - $m20*$m11\) | bc -l)
i21=$(echo -$invdet*\($m21*$m00 - $m20*$m01\) | bc -l)
i22=$(echo  $invdet*\($m11*$m00 - $m10*$m01\) | bc -l)

echo "xyz to rgb"
echo $i00 $i01 $i02
echo $i10 $i11 $i12
echo $i20 $i21 $i22

# this outputs gamma values for rgb:
read -r gr gg gb <<< $(iccdump -v3 -t rTRC -t gTRC -t bTRC $1 | grep gamma | awk '{print $NF}')
echo "gamma red, green, blue:"
echo $gr $gg $gb

name=custom

cat > include/colour/${name}.h<< EOF
static inline void colour_xyz_to_${name}(const float *const xyz, float *rgb)
{
  const float M[] =
  {
    $i00, $i01, $i02,
    $i10, $i11, $i12,
    $i20, $i21, $i22,
  };

  mat3_mulv(M, xyz, rgb);
  
  // apply tonecurve
  rgb[0] = powf(rgb[0], 1.f/$gr);
  rgb[1] = powf(rgb[1], 1.f/$gg);
  rgb[2] = powf(rgb[2], 1.f/$gb);
}

static inline void colour_${name}_to_xyz(const float *const rgb, float *xyz)
{
  const float M[] =
  {
    $m00, $m01, $m02,
    $m10, $m11, $m12,
    $m20, $m21, $m22,
  };

  float linear[3];
  // undo tonecurve
  linear[0] = powf(rgb[0], $gr);
  linear[1] = powf(rgb[1], $gg);
  linear[2] = powf(rgb[2], $gb);

  mat3_mulv(M, linear, xyz);
}

static inline void colour_${name}_print_info(FILE *f)
{
  fprintf(f, "$name from icc profile $1 adapted to whitepoint with chromaticity $x $y $z\n");
}
EOF
