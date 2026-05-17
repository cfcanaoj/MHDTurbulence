  
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
