#!/bin/bash

obj2geo=/home/jo/vcs/corona-13/tools/obj2geo
base="test_"
hair_scale=0.0005
target=/data/scenes/hairball

for i in $(seq 1 449)
do
  this=$(printf test_%06g.obj $i)
  next=$(printf test_%06g.obj $[$i+1])
  ${obj2geo} $this $next $hair_scale
  ind=$(printf %04g $i)
  for geo in *.geo
  do
    geobase=${geo%.geo}
    echo $geobase
    mv $geo $target/${geobase}_${ind}.geo
    cp $target/head.nra2 $target/anim_${ind}.nra2
    echo "2 ${geobase}_${ind}" >> $target/anim_${ind}.nra2
  done
done
