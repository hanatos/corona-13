#!/bin/bash

# lowres inner
rm voltest
CFLAGS="-DCLOUD_INNER -DVOLTEST_RES=0" make voltest
./voltest cloud-noise-inner-lo.vol -s

# highres outer
rm voltest
CFLAGS="-DCLOUD_OUTER -DVOLTEST_RES=1" make voltest
./voltest cloud-noise-outer.vol -s

# highres outer, cut in half 
rm voltest
CFLAGS="-DCLOUD_OUTER -DVOLTEST_RES=1 -DCLOUD_CUT" make voltest
./voltest cloud-noise-outer-cut.vol -s

# full combined cloud, for reference
rm voltest
CFLAGS="-DCLOUD_FULL -DVOLTEST_RES=1" make voltest
./voltest cloud-noise-both.vol -s

