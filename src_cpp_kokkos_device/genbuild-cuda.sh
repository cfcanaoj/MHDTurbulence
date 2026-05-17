#!/bin/bash

#####################
# Cuda
#####################

#cmake -B build_cuda -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON
#cmake -B build_cuda -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON
#cmake --build build_cuda -j

# --preset reads CMakePresets.json
cmake --preset cuda-ampere80
cmake --build --preset cuda-ampere80
cp  build/cuda-ampere80/Simulation.x ../exe
