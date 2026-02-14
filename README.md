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

To run the code, compile `Simulation.f90` on a GPU server. OpenACC version, `srcacc` is our defalut setup.
```bash
cd srcacc
make
```

Then `Simulation.x` is created in `../exe` directory.  
If you want to use a different version, go to the relevant directory shown below.

|directory|GPU/CPU|Language |Parallelization|
|:---|:---:|:---|
|`srcacc`    |GPU|Fortran|MPI OpenACC|
|`srcomp`    |GPU|Fortan |MPI OpenMP|
|`srcomp_cpp`|GPU|C++    |MPI OpenMP|
|`srccpu`    |CPU|Fortran|MPI OpenMP|

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
|`sj_xd.sh`|CPU|`srccpu`,|

---
## Data Output Specification

### 0. Performance Measurement Mode (No Intermediate Output)

For performance benchmarking, the code provides an option in `main.f90` to suppress intermediate outputs. In this case the initial and final snapshots are only damped.
```Fortran,
    logical,parameter::nooutput=.true.
```

The code supports three different output modes depending on the purpose of analysis

### 1. Quick check: Text output (ASCII)

If you want to quickly inspect the simulation results the code outputs text-format data. To employ this mode, edit `output.f90` as follows. 
```Fortran,
  logical,parameter:: binaryout= .false.
```
The data is damped as `ascdata/snap?????.csv` and format is "x y d vx vy p phi X". `gnuplot` is useful to quick check.
```bash
gnuplot
set view map
splot "ascdata/snap?????.csv" u 1:2:8 w pm3d
```

### 2. Full data visualization: XMF + binary (for VisIt / ParaView)
If you want to visualize the full 3D data, the code outputs binary data + XMF metadata. Check `output.f90` as follows. 
```Fortran,
  logical,parameter:: binaryout= .true.
```
The data is damped as `bindata/field?????.xmf`, `bindata/field?????.bin`, `bindata/grid1D.bin`, `bindata/grid2D.bin`, and `bindata/grid3D.bin`. Use `VisIt/ParaView` to check the data.

### 3. Detailed analysis: Analysis program
For more quantitative studies (spectra, statistics, etc.), use the analysis tools provided in the `analysis/` directory.
Check `output.f90` as follows. 
```Fortran,
  logical,parameter:: binaryout= .true.
```
The data is damped as `bindata/unf?????.dat`, `bindata/field?????.bin`, `bindata/grid1D.bin`, `bindata/grid2D.bin`, and `bindata/grid3D.bin`. 

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
