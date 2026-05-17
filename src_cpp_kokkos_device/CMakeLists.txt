cmake_minimum_required(VERSION 3.16)

project(MHDTurbulence LANGUAGES CXX)

add_executable(Simulation.x
    main.cpp
    mhd.cpp
    boundary.cpp
    output.cpp
    mpi_config.cpp
)

target_compile_features(Simulation.x PRIVATE cxx_std_17)

target_compile_options(Simulation.x PRIVATE
    -O3
    -march=native
    -Wall
)

find_package(MPI REQUIRED)
target_link_libraries(Simulation.x PRIVATE MPI::MPI_CXX)

find_package(OpenMP REQUIRED)
target_link_libraries(Simulation.x PRIVATE OpenMP::OpenMP_CXX)

include(FetchContent)
FetchContent_Declare(
  Kokkos
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.7.03.zip
  #URL https://github.com/kokkos/kokkos/archive/refs/tags/5.1.1.zip
)
FetchContent_MakeAvailable(Kokkos)

target_link_libraries(Simulation.x PRIVATE Kokkos::kokkos)
