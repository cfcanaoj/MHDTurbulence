#!/bin/bash

#####################
# OpenMP TARGET
#####################

#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_AMD_ZEN5=ON
#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_ARMV9_GRACE=ON
#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON 
#cmake --build build_openmp -j

# --preset reads CMakePresets.json
cmake --preset openmp-native
cmake --build --preset openmp-native
cp  build/openmp-native/Simulation.x ../exe
