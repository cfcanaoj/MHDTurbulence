# MHDTurbulence

This repository contains a 3D MHD solver (MPI + OpenACC/OpenMP target/OpenMP).
The current default test problem is the [**Kelvin–Helmholtz instability**](./KH.md).

<img width="320" height="240" alt="initial" src="https://github.com/user-attachments/assets/90d21b9e-e3f7-4606-b425-8d4d8715aab3" />
<img width="320" height="240" alt="final" src="https://github.com/user-attachments/assets/01ecbb59-c21b-4e38-8e6a-78d9e7433791" />

---

## How to copy the source code

After you login the server, copy the programs.
```bash
git clone git@github.com:cfcanaoj/MHDTurbulence MHDTurbulence
cd MHDTurbulence
```

---

## How to run

### Compile

To run the code, first compile source on a GPU server. The OpenACC version (`src_f90_acc_device`) is the default and recommended setup.
```bash
cd src_f90_acc_device
make
```

The executable `Simulation.x` is created in `../exe` directory.  
If you prefer another implementation, use the appropriate directory and codes listed below.

|directory|Language |GPU/CPU|Parallelization|
|:---|:---:|:---|:---|
|`src_f90_acc_device`    |Fortran|GPU|MPI OpenACC|
|`src_f90_omp_device`    |Fortran|GPU|MPI OpenMP Target|
|`src_f90_omp_host`      |Fortran|CPU|MPI OpenMP|
|`src_cpp_omp_device`    |C++    |GPU|MPI OpenMP|
|`src_cpp_kokkos_device` |C++    |GPU|MPI Kokkos|
|`src_cpp_kokkos_host`   |C++    |CPU|MPI Kokkos|

For stable simulations, please use the OpenACC (`src_f90_acc_device`) or CPU (`src_f90_omp_host`) versions.

### Execution
Copy the batch script and submit the job. The script for Slurm is prepared.
```bash
cd ../exe
cp ../misc/sj_g00.sh .
sbatch sj_g00.sh
```
Batch scripts depend on your parallelization scheme and environment.

|script|GPU/CPU|directory|
|:---|:---:|:---|
|`sj_g00.sh`|GPU|`src_f90_acc_device`, `src_f90_omp_device`, `src_cpp_omp_device`, `src_cpp_kokkos_device`|
|`sj_xd.sh`|CPU|`src_f90_omp_host`, `src_cpp_kokkos_host`|
