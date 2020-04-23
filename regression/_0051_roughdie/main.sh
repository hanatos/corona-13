#!/bin/bash

reflect=0
count=6
header="rough dielectric with phong lobe"

# first argument is roughness, then reflect or transmit (1 or 0)
echo "1.7 73 \#" | ../../tools/battle-test ../../shaders/libroughdielectric.so 0.10 $reflect $count > output.log

# create output, keep env variables
source ../makebattletest.sh
