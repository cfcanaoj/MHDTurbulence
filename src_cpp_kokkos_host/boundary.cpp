
/**
 * @file boundary.cpp
 * @brief 
 * @author Keigo Nitadori
 * @date 2025-08-21
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
#include <omp.h>

#include "config.hpp"
#include "mhd.hpp"
#include "boundary.hpp"

#include "mpi_config.hpp"

using namespace hydflux_mod;

namespace boundary_mod {
  // Default boundary configuration (same as the reference Fortran boundary.f90)
  int boundary_xin  = config::boundary_xin;
  int boundary_xout = config::boundary_xout;
  int boundary_yin  = config::boundary_yin;
  int boundary_yout = config::boundary_yout;
  int boundary_zin  = config::boundary_zin;
  int boundary_zout = config::boundary_zout;

  BoundaryArray<double> Bs,Br;


void AllocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace resolution_mod;

  int dev = omp_get_default_device();
  (void)dev;
  
  // buffer for send
  Bs.allocate(nprim ,ngh,ktot,jtot,itot);

  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<ngh; i++) {
	  Bs.Xs(n,k,j,i) = 0.0;
	  Bs.Xe(n,k,j,i) = 0.0;
  }
  
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<ngh; j++)
	for (int i=0; i<itot; i++) {
	  Bs.Ys(n,k,j,i) = 0.0;
	  Bs.Ye(n,k,j,i) = 0.0;
  }
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  Bs.Zs(n,k,j,i) = 0.0;
	  Bs.Ze(n,k,j,i) = 0.0;
  }

  // buffer for receive

  Br.allocate(nprim ,ngh,ktot,jtot,itot);

  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<ngh; i++) {
	  Br.Xs(n,k,j,i) = 0.0;
	  Br.Xe(n,k,j,i) = 0.0;
  }
  
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<ngh; j++)
	for (int i=0; i<itot; i++) {
	  Br.Ys(n,k,j,i) = 0.0;
	  Br.Ye(n,k,j,i) = 0.0;
  }
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  Br.Zs(n,k,j,i) = 0.0;
	  Br.Ze(n,k,j,i) = 0.0;
  }
};


void DeallocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace resolution_mod;
  using namespace hydflux_mod;

  Bs.deallocate();
  Br.deallocate();
}


void SendRecvBoundary(const BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace mpi_config_mod;
  using namespace resolution_mod;

  const int dev = omp_get_default_device();
  int rc;
  int nreq = 0;

  (void)dev, (void)rc;

  // ---- X direction ----
  if (ntiles[dir1] == 1) {
    const int bc_in  = boundary_xin;
    const int bc_out = boundary_xout;

    // Br.Xs : ghost at x-in  (left)
#if 0
    if (bc_in == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(n,k,j,i) = Bs.Xs(n,k,j,i);   // from x-out send buffer
    } else if (bc_in == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,ngh-1-i);
      // flip normal velocity v1
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ktot; k++)
        for (int j=0; j<jtot; j++)
          for (int i=0; i<ngh; i++)
            Br.Xs(nve1,k,j,i) = -Br.Xs(nve1,k,j,i);
    } else if (bc_in == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,0);
    }
#else
	Kokkos::parallel_for("BR_XS",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ngh, jtot, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (bc_in == periodicb) {
			for (int n=0; n<nprim; n++){
				Br.Xs(n,k,j,i) = Bs.Xs(n,k,j,i);   // from x-out send buffer
			}
		} else if (bc_in == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,ngh-1-i);
			}
			Br.Xs(nve1,k,j,i) = -Br.Xs(nve1,k,j,i);
		} else if (bc_in == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,0);
			}
		}
	});
#endif

    // Br.Xe : ghost at x-out (right)
#if 0
    if (bc_out == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(n,k,j,i) = Bs.Xe(n,k,j,i);   // from x-in send buffer
    } else if (bc_out == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1-i);
      // flip normal velocity v1
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ktot; k++)
        for (int j=0; j<jtot; j++)
          for (int i=0; i<ngh; i++)
            Br.Xe(nve1,k,j,i) = -Br.Xe(nve1,k,j,i);
    } else if (bc_out == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1);
    }
#else
	Kokkos::parallel_for("BR_XE",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ngh, jtot, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (bc_out== periodicb) {
			for (int n=0; n<nprim; n++){
				Br.Xe(n,k,j,i) = Bs.Xe(n,k,j,i);   // from x-in send buffer
			}
		} else if (bc_out== reflection) {
			for (int n=0; n<nprim; n++){
				Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1-i);
			}
			Br.Xe(nve1,k,j,i) = -Br.Xe(nve1,k,j,i);
		} else if (bc_out== outflow) {
			for (int n=0; n<nprim; n++){
				Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1);
			}
		}
	});
#endif
  } else {
    // MPI exchange where neighbors exist; apply reflection/outflow when neighbor is MPI_PROC_NULL
    // Ensure host copies of send buffers are up-to-date (Bs is produced on device)

    // MPI uses host pointers only
    double* h_Bs_Xe = Bs.Xe_data;
    double* h_Bs_Xs = Bs.Xs_data;
    double* h_Br_Xs = Br.Xs_data;
    double* h_Br_Xe = Br.Xe_data;
    if (n1m != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Xs, Br.size1, MPI_DOUBLE, n1m, 1100, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Xe, Bs.size1, MPI_DOUBLE, n1m, 1200, comm3d, &req[nreq++]);
    } else {
      // x-in physical boundary
#if 0
      if (boundary_xin == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<ngh; i++)
                Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,ngh-1-i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(nve1,k,j,i) = -Br.Xs(nve1,k,j,i);
      } else if (boundary_xin == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<ngh; i++)
                Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,0);
      }
#else
	Kokkos::parallel_for("BR_XS2",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ngh, jtot, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (boundary_xin == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,ngh-1-i);
			}
			Br.Xs(nve1,k,j,i) = -Br.Xs(nve1,k,j,i);
		} else if (boundary_xin == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,0);
			}
		}
	});
#endif
    }

    if (n1p != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Xe, Br.size1, MPI_DOUBLE, n1p, 1200, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Xs, Bs.size1, MPI_DOUBLE, n1p, 1100, comm3d, &req[nreq++]);
    } else {
      // x-out physical boundary
#if 0
      if (boundary_xout == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<ngh; i++)
                Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1-i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(nve1,k,j,i) = -Br.Xe(nve1,k,j,i);
      } else if (boundary_xout == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<ngh; i++)
                Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1);
      }
#else
	Kokkos::parallel_for("BR_XE2",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ngh, jtot, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (boundary_xout == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1-i);
			}
			Br.Xe(nve1,k,j,i) = -Br.Xe(nve1,k,j,i);
		} else if (boundary_xout == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1);
			}
		}
	});
#endif
    }
  }

  // ---- Y direction ----
  if (ntiles[dir2] == 1) {
    const int bc_in  = boundary_yin;
    const int bc_out = boundary_yout;

#if 0
    if (bc_in == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ys(n,k,j,i) = Bs.Ys(n,k,j,i);
    } else if (bc_in == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ys(n,k,j,i) = Bs.Ye(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ktot; k++)
        for (int j=0; j<ngh; j++)
          for (int i=0; i<itot; i++)
            Br.Ys(nve2,k,j,i) = -Br.Ys(nve2,k,j,i);
    } else if (bc_in == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ys(n,k,j,i) = Bs.Ye(n,k,0,i);
    }
#else
	Kokkos::parallel_for("BR_YS",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, ngh, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (bc_in == periodicb) {
			for (int n=0; n<nprim; n++){
				Br.Ys(n,k,j,i) = Bs.Ys(n,k,j,i);
			}
		} else if (bc_in == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Ys(n,k,j,i) = Bs.Ye(n,k,ngh-1-j,i);
			}
			Br.Ys(nve2,k,j,i) = -Br.Ys(nve2,k,j,i);
		} else if (bc_in == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Ys(n,k,j,i) = Bs.Ye(n,k,0,i);
			}
		}
	});
#endif

#if 0
    if (bc_out == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ye(n,k,j,i) = Bs.Ye(n,k,j,i);
    } else if (bc_out == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ktot; k++)
        for (int j=0; j<ngh; j++)
          for (int i=0; i<itot; i++)
            Br.Ye(nve2,k,j,i) = -Br.Ye(nve2,k,j,i);
    } else if (bc_out == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1,i);
    }
#else
	Kokkos::parallel_for("BR_YE",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, ngh, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (bc_out == periodicb) {
			for (int n=0; n<nprim; n++){
				Br.Ye(n,k,j,i) = Bs.Ye(n,k,j,i);
			}
		} else if (bc_out == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1-j,i);
			}
			Br.Ye(nve2,k,j,i) = -Br.Ye(nve2,k,j,i);
		} else if (bc_out == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1,i);
			}
		}
	});
#endif
  } else {
    // Ensure host copies of send buffers are up-to-date (Bs is produced on device)

    // MPI uses host pointers only
    double* h_Bs_Ye = Bs.Ye_data;
    double* h_Bs_Ys = Bs.Ys_data;
    double* h_Br_Ys = Br.Ys_data;
    double* h_Br_Ye = Br.Ye_data;
if (n2m != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Ys, Br.size2, MPI_DOUBLE, n2m, 2100, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Ye, Bs.size2, MPI_DOUBLE, n2m, 2200, comm3d, &req[nreq++]);
    } else {
#if 0
      if (boundary_yin == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<ngh; j++)
              for (int i=0; i<itot; i++)
                Br.Ys(n,k,j,i) = Bs.Ye(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ys(nve2,k,j,i) = -Br.Ys(nve2,k,j,i);
      } else if (boundary_yin == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<ngh; j++)
              for (int i=0; i<itot; i++)
                Br.Ys(n,k,j,i) = Bs.Ye(n,k,0,i);
      }
#else
	Kokkos::parallel_for("BR_YS2",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, ngh, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (boundary_yin == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Ys(n,k,j,i) = Bs.Ye(n,k,ngh-1-j,i);
			}
			Br.Ys(nve2,k,j,i) = -Br.Ys(nve2,k,j,i);
		} else if (boundary_yin == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Ys(n,k,j,i) = Bs.Ye(n,k,0,i);
			}
		}
	});
#endif
    }

    if (n2p != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Ye, Br.size2, MPI_DOUBLE, n2p, 2200, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Ys, Bs.size2, MPI_DOUBLE, n2p, 2100, comm3d, &req[nreq++]);
    } else {
#if 0
      if (boundary_yout == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<ngh; j++)
              for (int i=0; i<itot; i++)
                Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ktot; k++)
          for (int j=0; j<ngh; j++)
            for (int i=0; i<itot; i++)
              Br.Ye(nve2,k,j,i) = -Br.Ye(nve2,k,j,i);
      } else if (boundary_yout == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ktot; k++)
            for (int j=0; j<ngh; j++)
              for (int i=0; i<itot; i++)
                Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1,i);
      }
#else
	Kokkos::parallel_for("BR_YE2",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, ngh, ktot}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (boundary_yout == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1-j,i);
			}
			Br.Ye(nve2,k,j,i) = -Br.Ye(nve2,k,j,i);
		} else if (boundary_yout == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1,i);
			}
		}
	});
#endif
    }
  }

  // ---- Z direction ----
  if (ntiles[dir3] == 1) {
    const int bc_in  = boundary_zin;
    const int bc_out = boundary_zout;

#if 0
    if (bc_in == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Zs(n,k,j,i) = Bs.Zs(n,k,j,i);
    } else if (bc_in == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Zs(n,k,j,i) = Bs.Ze(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ngh; k++)
        for (int j=0; j<jtot; j++)
          for (int i=0; i<itot; i++)
            Br.Zs(nve3,k,j,i) = -Br.Zs(nve3,k,j,i);
    } else if (bc_in == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Zs(n,k,j,i) = Bs.Ze(n,0,j,i);
    }
#else
	Kokkos::parallel_for("BR_ZS",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, jtot, ngh}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (bc_in == periodicb) {
			for (int n=0; n<nprim; n++){
				Br.Zs(n,k,j,i) = Bs.Zs(n,k,j,i);
			}
		} else if (bc_in == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Zs(n,k,j,i) = Bs.Ze(n,ngh-1-k,j,i);
			}
			Br.Zs(nve3,k,j,i) = -Br.Zs(nve3,k,j,i);
		} else if (bc_in == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Zs(n,k,j,i) = Bs.Ze(n,0,j,i);
			}
		}
	});
#endif

#if 0
    if (bc_out == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Ze(n,k,j,i) = Bs.Ze(n,k,j,i);
    } else if (bc_out == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ngh; k++)
        for (int j=0; j<jtot; j++)
          for (int i=0; i<itot; i++)
            Br.Ze(nve3,k,j,i) = -Br.Ze(nve3,k,j,i);
    } else if (bc_out == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1,j,i);
    }
#else
	Kokkos::parallel_for("BR_ZE",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, jtot, ngh}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (bc_out== periodicb) {
			for (int n=0; n<nprim; n++){
				Br.Ze(n,k,j,i) = Bs.Ze(n,k,j,i);
			}
		} else if (bc_out== reflection) {
			for (int n=0; n<nprim; n++){
				Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1-k,j,i);
			}
			Br.Ze(nve3,k,j,i) = -Br.Ze(nve3,k,j,i);
		} else if (bc_out== outflow) {
			for (int n=0; n<nprim; n++){
				Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1,j,i);
			}
		}
	});
#endif
  } else {
    // Ensure host copies of send buffers are up-to-date (Bs is produced on device)

    // MPI uses host pointers only
    double* h_Bs_Ze = Bs.Ze_data;
    double* h_Bs_Zs = Bs.Zs_data;
    double* h_Br_Zs = Br.Zs_data;
    double* h_Br_Ze = Br.Ze_data;
if (n3m != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Zs, Br.size3, MPI_DOUBLE, n3m, 3100, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Ze, Bs.size3, MPI_DOUBLE, n3m, 3200, comm3d, &req[nreq++]);
    } else {
#if 0
      if (boundary_zin == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<itot; i++)
                Br.Zs(n,k,j,i) = Bs.Ze(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Zs(nve3,k,j,i) = -Br.Zs(nve3,k,j,i);
      } else if (boundary_zin == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<itot; i++)
                Br.Zs(n,k,j,i) = Bs.Ze(n,0,j,i);
      }
#else
	Kokkos::parallel_for("BR_ZS2",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, jtot, ngh}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (boundary_zin == reflection) {
			for (int n=0; n<nprim; n++){
				Br.Zs(n,k,j,i) = Bs.Ze(n,ngh-1-k,j,i);
			}
			Br.Zs(nve3,k,j,i) = -Br.Zs(nve3,k,j,i);
		} else if (boundary_zin == outflow) {
			for (int n=0; n<nprim; n++){
				Br.Zs(n,k,j,i) = Bs.Ze(n,0,j,i);
			}
		}
	});
#endif
    }

    if (n3p != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Ze, Br.size3, MPI_DOUBLE, n3p, 3200, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Zs, Bs.size3, MPI_DOUBLE, n3p, 3100, comm3d, &req[nreq++]);
    } else {
#if 0
      if (boundary_zout == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<itot; i++)
                Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ngh; k++)
          for (int j=0; j<jtot; j++)
            for (int i=0; i<itot; i++)
              Br.Ze(nve3,k,j,i) = -Br.Ze(nve3,k,j,i);
      } else if (boundary_zout == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=0; j<jtot; j++)
              for (int i=0; i<itot; i++)
                Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1,j,i);
      }
#else
	Kokkos::parallel_for("BR_ZE2",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {itot, jtot, ngh}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		if (boundary_zout== reflection) {
			for (int n=0; n<nprim; n++){
				Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1-k,j,i);
			}
			Br.Ze(nve3,k,j,i) = -Br.Ze(nve3,k,j,i);
		} else if (boundary_zout== outflow) {
			for (int n=0; n<nprim; n++){
				Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1,j,i);
			}
		}
	});
#endif
    }
  }

  
if (nreq != 0) MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);

#if 0
// After MPI completion, push received halo buffers to device before they are used.
// Only update faces that were actually received (avoid overwriting physical BC computed on device).
if (ntiles[dir1] != 1) {
  if (n1m != MPI_PROC_NULL) {
	#pragma omp target update to(Br.Xs_data[0:Br.size1])
	  }
  if (n1p != MPI_PROC_NULL) {
	#pragma omp target update to(Br.Xe_data[0:Br.size1])
	  }
}
if (ntiles[dir2] != 1) {
  if (n2m != MPI_PROC_NULL) {
	#pragma omp target update to(Br.Ys_data[0:Br.size2])
	  }
  if (n2p != MPI_PROC_NULL) {
	#pragma omp target update to(Br.Ye_data[0:Br.size2])
	  }
}
if (ntiles[dir3] != 1) {
  if (n3m != MPI_PROC_NULL) {
	#pragma omp target update to(Br.Zs_data[0:Br.size3])
	  }
  if (n3p != MPI_PROC_NULL) {
	#pragma omp target update to(Br.Ze_data[0:Br.size3])
	  }
}
#endif

nreq = 0;
};


void SetBoundaryCondition(FieldArray<double>& P,BoundaryArray<double>& Bs,BoundaryArray<double>& Br) {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;

  // x-direction
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xe|
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<ngh; i++) {
	  Bs.Xs(n,k,j,i) = P(n,k,j,ie-ngh+1+i);
	  Bs.Xe(n,k,j,i) = P(n,k,j,is      +i);
	  //Bs.Xs_data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+ie-ngh+1+i];
	  //Bs.Xe_data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+is      +i];
  }
#else
	Kokkos::parallel_for("PackX",
	Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {ngh, nprim, jtot, ktot}),
	KOKKOS_LAMBDA(const int i, const int n, const int j, const int k) {
		Bs.Xs(n,k,j,i) = P(n,k,j,ie-ngh+1+i);
		Bs.Xe(n,k,j,i) = P(n,k,j,is      +i);
	});
#endif
  
  // y-direction
  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<ngh;j++)
	for (int i=0; i<itot; i++) {
	  Bs.Ys(n,k,j,i) = P(n,k,je-ngh+1+j,i);
	  Bs.Ye(n,k,j,i) = P(n,k,js      +j,i);
	  //Bs.Ys_data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+je-ngh+1+j)*itot+i];
	  //Bs.Ye_data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+js      +j)*itot+i];
  }
#else
	Kokkos::parallel_for("PackY",
	Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {itot, nprim, ngh, ktot}),
	KOKKOS_LAMBDA(const int i, const int n, const int j, const int k) {
		Bs.Ys(n,k,j,i) = P(n,k,je-ngh+1+j,i);
		Bs.Ye(n,k,j,i) = P(n,k,js      +j,i);
	});
#endif
  
  // z-direction
  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  Bs.Zs(n,k,j,i) = P(n,ke-ngh+1+k,j,i);
	  Bs.Ze(n,k,j,i) = P(n,ks      +k,j,i);
	  //Bs.Zs_data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ke-ngh+1+k)*jtot+j)*itot+i];
	  //Bs.Ze_data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ks      +k)*jtot+j)*itot+i];
  }
#else
	Kokkos::parallel_for("PackZ",
	Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {itot, nprim, jtot, ngh}),
	KOKKOS_LAMBDA(const int i, const int n, const int j, const int k) {
		Bs.Zs(n,k,j,i) = P(n,ke-ngh+1+k,j,i);
		Bs.Ze(n,k,j,i) = P(n,ks      +k,j,i);
	});
#endif

  SendRecvBoundary(Bs,Br);

  
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xe|
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<ngh; i++) {
	  P(n,k,j,is-ngh+i) = Br.Xs(n,k,j,i);
	  P(n,k,j,ie+1  +i) = Br.Xe(n,k,j,i);
	  //P.data[((n*ktot+k)*jtot+j)*itot+is-ngh+i] = Xs.data[((n*ktot+k)*jtot+j)*ngh+i];
	  //P.data[((n*ktot+k)*jtot+j)*itot+ie+1  +i] = Xe.data[((n*ktot+k)*jtot+j)*ngh+i];
  }
#else
	Kokkos::parallel_for("UnpackX",
	Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {ngh, nprim, jtot, ktot}),
	KOKKOS_LAMBDA(const int i, const int n, const int j, const int k) {
		P.href(n,k,j,is-ngh+i) = Br.Xs(n,k,j,i);
		P.href(n,k,j,ie+1  +i) = Br.Xe(n,k,j,i);
	});
#endif

  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<ngh;j++)
	for (int i=0; i<itot; i++) {
	  P(n,k,js-ngh+j,i) = Br.Ys(n,k,j,i);
	  P(n,k,je+1  +j,i) = Br.Ye(n,k,j,i);
	  //P.data[((n*ktot+k)*jtot+js-ngh+j)*itot+i] = Ys.data[((n*ktot+k)*ngh+j)*itot+i];
	  //P.data[((n*ktot+k)*jtot+je+1  +j)*itot+i] = Ye.data[((n*ktot+k)*ngh+j)*itot+i];
  }
#else
	Kokkos::parallel_for("UnpackY",
	Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {itot, nprim, ngh, ktot}),
	KOKKOS_LAMBDA(const int i, const int n, const int j, const int k) {
		P.href(n,k,js-ngh+j,i) = Br.Ys(n,k,j,i);
		P.href(n,k,je+1  +j,i) = Br.Ye(n,k,j,i);
	});
#endif

  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  P(n,ks-ngh+k,j,i) = Br.Zs(n,k,j,i);
	  P(n,ke+1  +k,j,i) = Br.Ze(n,k,j,i);
	  //P.data[((n*ktot+ks-ngh+k)*jtot+j)*itot+i] = Zs.data[((n*ngh+k)*jtot+j)*itot+i];
	  //P.data[((n*ktot+ke+1  +k)*jtot+j)*itot+i] = Ze.data[((n*ngh+k)*jtot+j)*itot+i];
  }
#else
	Kokkos::parallel_for("UnpackZ",
	Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {itot, nprim, jtot, ngh}),
	KOKKOS_LAMBDA(const int i, const int n, const int j, const int k) {
		P.href(n,ks-ngh+k,j,i) = Br.Zs(n,k,j,i);
		P.href(n,ke+1  +k,j,i) = Br.Ze(n,k,j,i);
	});
#endif
};

};// end namespace
