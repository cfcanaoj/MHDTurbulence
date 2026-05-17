#!/bin/bash


#####################
# OpenMP
#####################

#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_AMD_ZEN5=ON
#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_ARMV9_GRACE=ON
#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON 
#cmake --build build_openmp -j

# --preset reads CMakePresets.json
cmake --preset openmp-native
cmake --build --preset openmp-native

#####################
# Cuda
#####################

#cmake -B build_cuda -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON
#cmake -B build_cuda -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON
#cmake --build build_cuda -j

# --preset reads CMakePresets.json
cmake --preset cuda-ampere80
cmake --build --preset cuda-ampere80
