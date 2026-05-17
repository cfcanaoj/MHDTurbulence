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

Boundary conditions are applied in `config.f90`. In the default, reflection boundary is taken for y-direction and periodic boundart is taken for other directions.
```Fortran
  integer,parameter:: periodicb=1,reflection=2,outflow=3
  integer,parameter:: boundary_xin=periodicb , boundary_xout=periodicb
  integer,parameter:: boundary_yin=reflection, boundary_yout=reflection
  integer,parameter:: boundary_zin=periodicb , boundary_zout=periodicb
```
