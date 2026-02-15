# MHDTurbulence

This repository contains a 3D MHD solver (MPI + OpenACC/OpenMP target/OpenMP).
The current default test problem is the [**Kelvin–Helmholtz instability**](./KH.md).

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

To run the code, first compile source on a GPU server. The OpenACC version (`srcacc`) is the defalut and recommended setup.
```bash
cd srcacc
make
```

The executable `Simulation.x` is created in `../exe` directory.  
If you prefer another implementation, use the appropriate directory and codes listed below.

|directory|Language |GPU/CPU|Parallelization|
|:---|:---:|:---|:---|
|`srcacc`    |Fortran|GPU|MPI OpenACC|
|`srcomp`    |Fortan |GPU|MPI OpenMP|
|`srcomp_cpp`|C++    |GPU|MPI OpenMP|
|`srccpu`    |Fortran|CPU|MPI OpenMP|

The OpenMP version (`srcomp`) is currently not fully functional and may not run correctly.
For stable simulations, please use the OpenACC (`srcacc`) or CPU
(`srccpu`) versions.

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
|`sj_g00.sh`|GPU|`srcacc`, `srcomp`, `srcomp_cpp`|
|`sj_xd.sh`|CPU|`srccpu`|

---
## Data Output Specification

The code supports multiple output modes depending on your purpose:
benchmarking, quick inspection, visualization, or detailed analysis.

### 0. Performance Measurement Mode (No Intermediate Output)

For performance benchmarking, the code provides an option in `config.f90` to suppress intermediate outputs. In this case the initial and final snapshots are only damped.
```Fortran,
logical,parameter:: benchmarkmode = .true. !! If true, only initial and final outputs are damped.
```

When enabled:

-   Intermediate snapshots are not written
-   I/O overhead is minimized
-   Only initial and final data are output

This mode is recommended for:

-   Strong/weak scaling tests
-   GPU/CPU performance comparisons
-   Pure solver benchmarking

#### Performance Information
During execution, the following performance information is printed to the standard output:
```bash
sim time [s]: 1.145357e+03
time/count/cell : 2.262434e-09
```
- sim time [s] : total wall-clock time of the main loop of simulation
- time/count/cell : wall-clock time per cell per time step
These values are useful for benchmarking and performance comparison.

#### Benchmark Results
The typical performance in representative environments is shown below.
- Grid size: number of cells in each direction
- Physicsl time (t): physical end time of the simulation
- Wall time: total elapsed wall-clock time
- time/cell/step: wall-clock time per cell per time step

|Code|Grid size x Physical time|Wall time [s]|time/cell/step [s]|Environment|
|:---|:---:|---:|---:|:---|
|`srcacc`    |150^3 x 15|408.18|8.06e-10|CfCA GPU server, A100 4 GPU|
|`srcomp_cpp`|150^3 x 15|1145  |22.6e-10|CfCA GPU server, A100 4 GPU|
|`srccpu`    |156^3 x 15|      |        |CfCA XD2000, Xeon Max 1 node|


### 1. Quick check: Text output (ASCII)

For quick inspection and debugging, ASCII output can be enabled in
`config.f90`:
```Fortran,
logical,parameter:: asciiout = .true. !! Ascii-files are additionaly damped.
```
The data is damped as `ascdata/snap###-?????.csv`. Here ### is the mpi-process and ????? is the number of the snapshots. The format is "x y d vx vy p phi X". `gnuplot` is useful to quick check.
```bash
gnuplot
set view map
splot "ascdata/snap###-?????.csv" u 1:2:8 w pm3d
```
You can compare your data with [the sample](./sampledata). Since we use random perturbation. The pattern of the turbulence is not necessarily match the sample.

<img width="320" height="240" alt="initial" src="https://github.com/user-attachments/assets/90d21b9e-e3f7-4606-b425-8d4d8715aab3" />
<img width="320" height="240" alt="final" src="https://github.com/user-attachments/assets/01ecbb59-c21b-4e38-8e6a-78d9e7433791" />


### 2. Full data visualization: XMF + binary (for VisIt / ParaView)
For full 3D visualization, the code outputs binary data together with
XMF metadata files.

Generated files:

    bindata/field?????.xmf
    bindata/field?????.bin
    bindata/grid1D.bin
    bindata/grid2D.bin
    bindata/grid3D.bin

Open the `.xmf` file directly in:

-   VisIt
-   ParaView

### 3. Detailed analysis: Analysis program
For more quantitative studies (spectra, statistics, etc.), use the analysis tools provided in the `analysis/` directory. The data is damped as `bindata/unf?????.dat`, `bindata/field?????.bin`, `bindata/grid1D.bin`, `bindata/grid2D.bin`, and `bindata/grid3D.bin`. 

#### Build the analysis tool

```bash
cd ../analysis
ln -s ../exe/bindata .
make
```

Count time snapshots:

```bash
./CountBindata.sh
cat control.dat
```

Run analysis:

```bash
sbatch sj_g00_ana.sh
```

Outputs are saved in `output/`.

#### 2D plots and animation

Copy `bindata/` and `output/` to the analysis server and:

```bash
make 2Dsnaps
make movie
```

Images are saved in `figures/`, movies in `movie/anivor`.

#### Spectrum

```bash
make spectrum
```

---
