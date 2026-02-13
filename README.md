# MHDTurbulence:

This repository contains a 3D MHD solver (MPI + OpenACC/OpenMP/CPU versions) and example setups.
The current default problem is **Kelvin–Helmholtz instability**.

---

## How to copy the source code

After you login the server, `g00.cfca.nao.ac.jp`, follow the instruction.

```bash
cd /gwork0/<username>
git clone git@github.com:cfcanaoj/MHDTurbulence MHDTurbulence
cd MHDTurbulence
```

---

## How to run

### Compile

To run the code, compile `Simulation.f90` on a GPU server.

```bash
cd srcacc
make
```

Then `Simulation.x` is created in `../exe` directory (OpenACC version).  
If you want to use a different version, go to the relevant directory shown below.

|directory|GPU/CPU|Parallelization|
|:---|:---:|:---|
|`srcacc`|GPU|MPI OpenACC|
|`srcomp`|GPU|MPI OpenMP|
|`srccpu`|CPU|MPI|

### Run

```bash
cd ../exe
cp ../misc/sj_g00.sh .
sbatch sj_g00.sh
```

The simulation data are saved in `bindata/`.

Batch scripts depend on your parallelization scheme and environment.

|script|GPU/CPU|Parallelization|
|:---|:---:|:---|
|`sj_g00.sh`|GPU|MPI OpenACC|
|`sj_g00.sh`|GPU|MPI OpenMP|
|`sj_xd.sh`|CPU|MPI|

---

## Problem setup: Kelvin–Helmholtz instability

The initial condition is defined in `GenerateProblem` in `main.f90`.

### Hydrodynamic variables

Uniform density and pressure with $$\rho = 1$$ and $$p = 1$$.
A shear flow in the **x-direction** (two shear layers centered at $$y=\pm 0.5$$):

- $$\Delta v = 2.0$$
- shear-layer width $$w = 0.05$$

$$
v_x(y)=\frac{\Delta v}{2}\left[\tanh\left(\frac{y+0.5}{w}\right)-\tanh\left(\frac{y-0.5}{w}\right)-1\right]
$$

A small perturbation in the **y-direction** to seed KH roll-up (sinusoidal in x, localized around the layers):

- amplitude $$10^{-3}$$
- localization width $$\sigma = 0.2$$

$$
v_y(x,y)=10^{-3}\sin(2\pi x)\left[e^{-(y+0.5)^2/\sigma^2}+e^{-(y-0.5)^2/\sigma^2}\right]
$$
$$
v_z = 0
$$

Additionally, a small **random perturbation** is added to $$v_x$$ (1% of $$\Delta v$$):

- $$\mathrm{rrv}=10^{-2}$$
- $$v_x \leftarrow v_x + \Delta v\,\mathrm{rrv}\,(r-0.5)\) where \(r\in[0,1)$$

### Magnetic field

In the current KH setup, the magnetic field is initialized to zero:

$$
\mathbf{B} = (0,0,0)
$$

### Tracer / composition field

`Xcomp(1,...)` is initialized as a smooth step-like tracer associated with the two layers, useful for visualizing mixing.

### Boundary conditions

Boundary conditions are applied via `BoundaryCondition` (see `boundarymod`).  
(The domain extent and grid size are configured in `basicmod` and printed at runtime.)

---

## Analysis

### Build the analysis tool

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

### 2D plots and animation

Copy `bindata/` and `output/` to the analysis server and:

```bash
make 2Dsnaps
make movie
```

Images are saved in `figures/`, movies in `movie/anivor`.

### Spectrum

```bash
make spectrum
```

---

# Results
Typical results are shown [here](https://www.notion.so/Turbulence-Studies-numerical-experiments-e4836ad642684f8f992d54a1f7e22635).
