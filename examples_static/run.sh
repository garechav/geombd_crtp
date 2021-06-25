#!/bin/bash
echo "Running forward dynamics examples"

cd build/examples
./example_FD_lwr
./example_FD_nao
./example_FD_hrp2
./example_FD_atlas

echo "Running forward dynamics differentiation examples"

./example_D_FD_lwr
./example_D_FD_nao
./example_D_FD_hrp2
./example_D_FD_atlas