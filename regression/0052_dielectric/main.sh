#!/bin/bash

reflect=1
count=4
header="rough dielectric reflection ggx/visible normals/smith"

# first argument is roughness, then reflect or transmit (1 or 0)
echo "1.7 73 \#" | ../../tools/battle-test ../../shaders/libdielectric.so 0.4 $reflect $count > output.log

# create output, keep env variables
source ../makebattletest.sh
