# MHDTurbulence:

This repository contains a 3D MHD solver (MPI + OpenACC/OpenMP/CPU versions) and example setups.
The current default problem is [**Kelvin–Helmholtz instability**](./KH.md).

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

To run the code, first compile the source codes in a GPU server. OpenACC version, `srcacc` is our defalut setup.
```bash
cd srcacc
make
```

Then `Simulation.x` is created in `../exe` directory.  
If you want to use a different version, go to the relevant directory shown below.

|directory|Language |GPU/CPU|Parallelization|
|:---|:---:|:---|:---|
|`srcacc`    |Fortran|GPU|MPI OpenACC|
|`srcomp`    |Fortan |GPU|MPI OpenMP|
|`srcomp_cpp`|C++    |GPU|MPI OpenMP|
|`srccpu`    |Fortran|CPU|MPI OpenMP|

The OpenMP version (`srcomp`) is currently not fully functional and may not run correctly.
Please use the OpenACC (GPU) or CPU version instead.

### Run
Copy the batch script and run the code. The script for Slurm is prepared.
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

### 0. Performance Measurement Mode (No Intermediate Output)

For performance benchmarking, the code provides an option in `config.f90` to suppress intermediate outputs. In this case the initial and final snapshots are only damped.
```Fortran,
logical,parameter:: benchmarkmode = .true. !! If true, only initial and final outputs are damped.
```

The code supports three different output modes depending on the purpose of analysis

### 1. Quick check: Text output (ASCII)

If you want to quickly inspect the simulation results the code outputs text-format data. To employ this mode, edit `config.f90` as follows. 
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
If you want to visualize the full 3D data, the code outputs binary data + XMF metadata. The data is damped as `bindata/field?????.xmf`, `bindata/field?????.bin`, `bindata/grid1D.bin`, `bindata/grid2D.bin`, and `bindata/grid3D.bin`. Use `VisIt/ParaView` to check the data.

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
