#!/bin/bash
echo "Cinvestav 2021"

echo "---------------------------------"
echo "Running forward dynamics examples"
echo "---------------------------------"

cd examples_static/build/examples
./example_FD lwr.urdf
./example_FD nao_inertial_python.urdf
./example_FD HRP2.urdf
./example_FD atlas.urdf

echo "-------------------------------------------------"
echo "Running forward dynamics differentiation examples"
echo "-------------------------------------------------"

./example_D_FD lwr.urdf
./example_D_FD nao_inertial_python.urdf
./example_D_FD HRP2.urdf
./example_D_FD atlas.urdf
