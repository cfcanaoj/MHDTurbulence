## Problem setup: Kelvin-Helmholtz instability

The initial condition is defined in `GenerateProblem` in `main.f90`.

### Hydrodynamic variables

Uniform density and pressure are set as

$$
\rho = 1,
\qquad
p = 1.
$$

A shear flow is imposed in the $x$-direction with two shear layers centered at $y=\pm 0.5$.
The main parameters are:

- $\Delta v = 2.0$
- shear-layer width $w = 0.05$

The initial $x$-velocity is

$$
v_x(y)=\frac{\Delta v}{2}\left[
\tanh\left(\frac{y+0.5}{w}\right)
-\tanh\left(\frac{y-0.5}{w}\right)
-1
\right].
$$

A small perturbation is added in the $y$-direction to seed KH roll-up.
It is sinusoidal in $x$ and localized around the two shear layers.
The parameters are:

- amplitude $10^{-2}$
- localization width $\sigma = 0.2$

The perturbation is

$$
v_y(x,y)=10^{-2}\sin(2\pi x)\left[
e^{-(y+0.5)^2/\sigma^2}
+
e^{-(y-0.5)^2/\sigma^2}
\right].
$$

The $z$-velocity is initially

$$
v_z = 0.
$$

Additionally, the code contains a random perturbation term for $v_x$, but in the current setup it is disabled:

- $\mathrm{rrv}=0.0\times 10^{-2}$

If enabled, the perturbation is

$$
v_x \leftarrow v_x + \Delta v\,\mathrm{rrv}\,(r-0.5),
\qquad r\in[0,1).
$$

### Magnetic field

In the current KH setup, the magnetic field is initialized to zero:

$$
\mathbf{B} = (0,0,0).
$$

### Tracer / composition field

`Xcomp(1,...)` is initialized as a smooth step-like tracer associated with the two layers, useful for visualizing mixing:

$$
X_{\mathrm{comp}} =
\frac{1}{2}\left[
\tanh\left(\frac{y+0.5}{w}\right)
- \tanh\left(\frac{y-0.5}{w}\right)
\right].
$$

### Boundary conditions

Boundary conditions are applied in `config.f90`. In the default setup, reflection boundary conditions are used in the $y$-direction and periodic boundary conditions are used in the other directions.

```Fortran
  integer,parameter:: periodicb=1,reflection=2,outflow=3
  integer,parameter:: boundary_xin=periodicb , boundary_xout=periodicb
  integer,parameter:: boundary_yin=reflection, boundary_yout=reflection
  integer,parameter:: boundary_zin=periodicb , boundary_zout=periodicb
```

### Real-time analysis

`RealTimeAnalysis` in `src_f90_acc_device/main.f90` evaluates bulk diagnostics during the run and writes them to `t-prof.csv`.

For each cell, the cell volume is computed as

$$
\Delta V_{ijk} =
\left(x1a_{i+1}-x1a_i\right)
\left(x2a_{j+1}-x2a_j\right)
\left(x3a_{k+1}-x3a_k\right).
$$

The routine then performs MPI reductions for

$$
V = \sum_{ijk}\Delta V_{ijk},
$$

$$
M = \sum_{ijk} X_{\mathrm{comp},ijk}\left(1-X_{\mathrm{comp},ijk}\right)\Delta V_{ijk},
$$

$$
S_{v_y} = \sum_{ijk} v_{y,ijk}^2\,\Delta V_{ijk}.
$$

The reported diagnostics are

$$
\mathrm{mix} = \frac{M}{V},
$$

and

$$
\left< v_y^2 \right>^{1/2} = \sqrt{\frac{S_{v_y}}{V}}.
$$

Here, `mix` is a simple mixing indicator based on the tracer field. It is close to zero in unmixed regions where $X_{\mathrm{comp}}$ is near 0 or 1, and becomes larger as the two layers mix.

The file `t-prof.csv` contains four columns:

1. `time`
2. `mix`
3. `sqrt(<v_y^2>)`
4. `A*exp(Gamma*time)`

The fourth column is a reference exponential growth curve:

The growth rate is $\Gamma$.

$$
A\exp(\Gamma t)
$$

with the hard-coded parameters

- $A = 1.2\times 10^{-3}$
- $\Gamma = 1.49$

These values are used as a practical comparison metric for KH growth in this setup. The source comment notes that this is not a general analytic growth rate for finite-width shear layers, but an empirical reference value for this problem.
