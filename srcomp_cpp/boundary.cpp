
/**
 * @file boundary.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
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

#pragma omp declare target
  BoundaryArray<double> Bs,Br;
#pragma omp end declare target 

auto assocb = [&](void* host_ptr, size_t bytes, int dev) {
    void* dptr = omp_get_mapped_ptr(host_ptr, dev);
    if (!dptr) {
        std::fprintf(stderr, "mapped_ptr is NULL for %p\n", host_ptr);
    }
    int r = omp_target_associate_ptr(host_ptr, dptr, bytes, /*device_offset=*/0, dev);
    if (r != 0) {
        std::fprintf(stderr, "omp_target_associate_ptr failed (%d)\n", r);
    }
};

void AllocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace resolution_mod;

  int dev = omp_get_default_device();
  
  // buffer for send
  Bs.allocate(nprim ,ngh,ktot,jtot,itot);

#pragma omp target enter data map (alloc: Bs.Xs_data[0:Bs.size1], Bs.Xe_data[0: Bs.size1])
#pragma omp target enter data map (alloc: Bs.Ys_data[0:Bs.size2], Bs.Ye_data[0: Bs.size2])
#pragma omp target enter data map (alloc: Bs.Zs_data[0:Bs.size3], Bs.Ze_data[0: Bs.size3])
#pragma omp target update to (Bs.n1, Bs.n2, Bs.n3, Bs.ng, Bs.nv, Bs.size1, Bs.size2, Bs.size3)
  /*
  assocb(Bs.Xs_data, sizeof(double)*Bs.size1,dev);
  assocb(Bs.Xe_data, sizeof(double)*Bs.size1,dev);
  assocb(Bs.Ys_data, sizeof(double)*Bs.size2,dev);
  assocb(Bs.Ye_data, sizeof(double)*Bs.size2,dev);
  assocb(Bs.Zs_data, sizeof(double)*Bs.size3,dev);
  assocb(Bs.Ze_data, sizeof(double)*Bs.size3,dev);
  */
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
#pragma omp target update to (Bs.Xs_data[0:Bs.size1],Bs.Xe_data[0:Bs.size1])
#pragma omp target update to (Bs.Ys_data[0:Bs.size2],Bs.Ye_data[0:Bs.size2])
#pragma omp target update to (Bs.Zs_data[0:Bs.size3],Bs.Ze_data[0:Bs.size3])

  // buffer for receive

  Br.allocate(nprim ,ngh,ktot,jtot,itot);

#pragma omp target enter data map (alloc: Br.Xs_data[0:Br.size1], Br.Xe_data[0: Br.size1])
#pragma omp target enter data map (alloc: Br.Ys_data[0:Br.size2], Br.Ye_data[0: Br.size2])
#pragma omp target enter data map (alloc: Br.Zs_data[0:Br.size3], Br.Ze_data[0: Br.size3])
#pragma omp target update to (Br.n1, Br.n2, Br.n3, Br.ng, Br.nv, Br.size1, Br.size2, Br.size3)
  /*
  assocb(Br.Xs_data, sizeof(double)*Br.size1,dev);
  assocb(Br.Xe_data, sizeof(double)*Br.size1,dev);
  assocb(Br.Ys_data, sizeof(double)*Br.size2,dev);
  assocb(Br.Ye_data, sizeof(double)*Br.size2,dev);
  assocb(Br.Zs_data, sizeof(double)*Br.size3,dev);
  assocb(Br.Ze_data, sizeof(double)*Br.size3,dev);
  */
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

#pragma omp target update to (Br.Xs_data[0:Br.size1],Br.Xe_data[0:Br.size1])
#pragma omp target update to (Br.Ys_data[0:Br.size2],Br.Ye_data[0:Br.size2])
#pragma omp target update to (Br.Zs_data[0:Br.size3],Br.Ze_data[0:Br.size3])


};


void DeallocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace resolution_mod;
  using namespace hydflux_mod;


#pragma omp target exit data map (delete: Bs.Xs_data[0:Bs.size1], Bs.Xe_data[0: Bs.size1])
#pragma omp target exit data map (delete: Bs.Ys_data[0:Bs.size2], Bs.Ye_data[0: Bs.size2])
#pragma omp target exit data map (delete: Bs.Zs_data[0:Bs.size3], Bs.Ze_data[0: Bs.size3])
  
#pragma omp target exit data map (delete: Br.Xs_data[0:Br.size1], Br.Xe_data[0: Br.size1])
#pragma omp target exit data map (delete: Br.Ys_data[0:Br.size2], Br.Ye_data[0: Br.size2])
#pragma omp target exit data map (delete: Br.Zs_data[0:Br.size3], Br.Ze_data[0: Br.size3])
}


void SendRecvBoundary(const BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace mpi_config_mod;
  using namespace resolution_mod;

	  const int dev  = omp_get_default_device();
	  const int host = omp_get_initial_device();

	  // If p is a mapped host pointer, return the device mapping.
	  // If p is already a device pointer (or not mapped), fall back to p.
	  auto devptr_ro = [&](const void* p)->const void*{
	    const void* dp = omp_get_mapped_ptr(const_cast<void*>(p), dev);
	    return dp ? dp : p;
	  };
	  auto devptr_rw = [&](void* p)->void*{
	    void* dp = omp_get_mapped_ptr(p, dev);
	    return dp ? dp : p;
	  };

	  // Dedicated host staging buffers (allocated once, reused). This guarantees MPI never
	  // touches device pointers even if the BoundaryArray storage changes in the future.
	  static std::vector<double> hsend_Xs, hsend_Xe, hrecv_Xs, hrecv_Xe;
	  static std::vector<double> hsend_Ys, hsend_Ye, hrecv_Ys, hrecv_Ye;
	  static std::vector<double> hsend_Zs, hsend_Ze, hrecv_Zs, hrecv_Ze;
  int rc;
  int nreq = 0;

  // ---- X direction ----
  if (ntiles[dir1] == 1) {
    const int bc_in  = boundary_xin;
    const int bc_out = boundary_xout;

    // Br.Xs : ghost at x-in  (left)
    if (bc_in == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(n,k,j,i) = Bs.Xs(n,k,j,i);   // from x-out send buffer
    } else if (bc_in == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,ngh-1-i);
      // flip normal velocity v1
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=ks; k<=ke; k++)
        for (int j=js; j<=je; j++)
          for (int i=0; i<ngh; i++)
            Br.Xs(nve1,k,j,i) = -Br.Xs(nve1,k,j,i);
    } else if (bc_in == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,0);
    }

    // Br.Xe : ghost at x-out (right)
    if (bc_out == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(n,k,j,i) = Bs.Xe(n,k,j,i);   // from x-in send buffer
    } else if (bc_out == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1-i);
      // flip normal velocity v1
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=ks; k<=ke; k++)
        for (int j=js; j<=je; j++)
          for (int i=0; i<ngh; i++)
            Br.Xe(nve1,k,j,i) = -Br.Xe(nve1,k,j,i);
    } else if (bc_out == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1);
    }
  } else {
    // MPI exchange where neighbors exist; apply reflection/outflow when neighbor is MPI_PROC_NULL
    // IMPORTANT:
    //   Do NOT pass device pointers to MPI unless you are 100% sure the MPI stack is
    //   CUDA-aware and correctly configured for your OpenMP offload runtime.
    //   We therefore stage through host buffers:
    //     device -> host (send buffers), MPI on host, host -> device (recv buffers).
    // This avoids cuMemGetAddressRange / CUDA-aware detection paths entirely.

	    // Stage device send buffers -> host staging buffers.
	    hsend_Xe.resize(Bs.size1);
	    hsend_Xs.resize(Bs.size1);
	    hrecv_Xs.resize(Br.size1);
	    hrecv_Xe.resize(Br.size1);
	    omp_target_memcpy(hsend_Xe.data(), devptr_ro(Bs.Xe_data), sizeof(double)*Bs.size1, 0, 0, host, dev);
	    omp_target_memcpy(hsend_Xs.data(), devptr_ro(Bs.Xs_data), sizeof(double)*Bs.size1, 0, 0, host, dev);

	    double* h_Bs_Xe = hsend_Xe.data();
	    double* h_Bs_Xs = hsend_Xs.data();
	    double* h_Br_Xs = hrecv_Xs.data();
	    double* h_Br_Xe = hrecv_Xe.data();

    if (n1m != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Xs, Br.size1, MPI_DOUBLE, n1m, 1100, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Xe, Bs.size1, MPI_DOUBLE, n1m, 1200, comm3d, &req[nreq++]);
    } else {
      // x-in physical boundary
      if (boundary_xin == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=js; j<=je; j++)
              for (int i=0; i<ngh; i++)
                Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,ngh-1-i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xs(nve1,k,j,i) = -Br.Xs(nve1,k,j,i);
      } else if (boundary_xin == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=js; j<=je; j++)
              for (int i=0; i<ngh; i++)
                Br.Xs(n,k,j,i) = Bs.Xe(n,k,j,0);
      }
    }

    if (n1p != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Xe, Br.size1, MPI_DOUBLE, n1p, 1200, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Xs, Bs.size1, MPI_DOUBLE, n1p, 1100, comm3d, &req[nreq++]);
    } else {
      // x-out physical boundary
      if (boundary_xout == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=js; j<=je; j++)
              for (int i=0; i<ngh; i++)
                Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1-i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=ks; k<=ke; k++)
          for (int j=js; j<=je; j++)
            for (int i=0; i<ngh; i++)
              Br.Xe(nve1,k,j,i) = -Br.Xe(nve1,k,j,i);
      } else if (boundary_xout == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=js; j<=je; j++)
              for (int i=0; i<ngh; i++)
                Br.Xe(n,k,j,i) = Bs.Xs(n,k,j,ngh-1);
      }
    }
  }

  // ---- Y direction ----
  if (ntiles[dir2] == 1) {
    const int bc_in  = boundary_yin;
    const int bc_out = boundary_yout;

    if (bc_in == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ys(n,k,j,i) = Bs.Ys(n,k,j,i);
    } else if (bc_in == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ys(n,k,j,i) = Bs.Ye(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=ks; k<=ke; k++)
        for (int j=0; j<ngh; j++)
          for (int i=is; i<=ie; i++)
            Br.Ys(nve2,k,j,i) = -Br.Ys(nve2,k,j,i);
    } else if (bc_in == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ys(n,k,j,i) = Bs.Ye(n,k,0,i);
    }

    if (bc_out == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ye(n,k,j,i) = Bs.Ye(n,k,j,i);
    } else if (bc_out == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=ks; k<=ke; k++)
        for (int j=0; j<ngh; j++)
          for (int i=is; i<=ie; i++)
            Br.Ye(nve2,k,j,i) = -Br.Ye(nve2,k,j,i);
    } else if (bc_out == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1,i);
    }
  } else {
	    // Host staging for MPI (see X direction block for rationale)
	    hsend_Ye.resize(Bs.size2);
	    hsend_Ys.resize(Bs.size2);
	    hrecv_Ys.resize(Br.size2);
	    hrecv_Ye.resize(Br.size2);
	    omp_target_memcpy(hsend_Ye.data(), devptr_ro(Bs.Ye_data), sizeof(double)*Bs.size2, 0, 0, host, dev);
	    omp_target_memcpy(hsend_Ys.data(), devptr_ro(Bs.Ys_data), sizeof(double)*Bs.size2, 0, 0, host, dev);

	    double* h_Bs_Ye = hsend_Ye.data();
	    double* h_Bs_Ys = hsend_Ys.data();
	    double* h_Br_Ys = hrecv_Ys.data();
	    double* h_Br_Ye = hrecv_Ye.data();

    if (n2m != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Ys, Br.size2, MPI_DOUBLE, n2m, 2100, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Ye, Bs.size2, MPI_DOUBLE, n2m, 2200, comm3d, &req[nreq++]);
    } else {
      if (boundary_yin == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=0; j<ngh; j++)
              for (int i=is; i<=ie; i++)
                Br.Ys(n,k,j,i) = Bs.Ye(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ys(nve2,k,j,i) = -Br.Ys(nve2,k,j,i);
      } else if (boundary_yin == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=0; j<ngh; j++)
              for (int i=is; i<=ie; i++)
                Br.Ys(n,k,j,i) = Bs.Ye(n,k,0,i);
      }
    }

    if (n2p != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Ye, Br.size2, MPI_DOUBLE, n2p, 2200, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Ys, Bs.size2, MPI_DOUBLE, n2p, 2100, comm3d, &req[nreq++]);
    } else {
      if (boundary_yout == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=0; j<ngh; j++)
              for (int i=is; i<=ie; i++)
                Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1-j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=ks; k<=ke; k++)
          for (int j=0; j<ngh; j++)
            for (int i=is; i<=ie; i++)
              Br.Ye(nve2,k,j,i) = -Br.Ye(nve2,k,j,i);
      } else if (boundary_yout == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=ks; k<=ke; k++)
            for (int j=0; j<ngh; j++)
              for (int i=is; i<=ie; i++)
                Br.Ye(n,k,j,i) = Bs.Ys(n,k,ngh-1,i);
      }
    }
  }

  // ---- Z direction ----
  if (ntiles[dir3] == 1) {
    const int bc_in  = boundary_zin;
    const int bc_out = boundary_zout;

    if (bc_in == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Zs(n,k,j,i) = Bs.Zs(n,k,j,i);
    } else if (bc_in == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Zs(n,k,j,i) = Bs.Ze(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ngh; k++)
        for (int j=js; j<=je; j++)
          for (int i=is; i<=ie; i++)
            Br.Zs(nve3,k,j,i) = -Br.Zs(nve3,k,j,i);
    } else if (bc_in == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Zs(n,k,j,i) = Bs.Ze(n,0,j,i);
    }

    if (bc_out == periodicb) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Ze(n,k,j,i) = Bs.Ze(n,k,j,i);
    } else if (bc_out == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
      for (int k=0; k<ngh; k++)
        for (int j=js; j<=je; j++)
          for (int i=is; i<=ie; i++)
            Br.Ze(nve3,k,j,i) = -Br.Ze(nve3,k,j,i);
    } else if (bc_out == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
      for (int n=0; n<nprim; n++)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1,j,i);
    }
  } else {
	    // Host staging for MPI (see X direction block for rationale)
	    hsend_Ze.resize(Bs.size3);
	    hsend_Zs.resize(Bs.size3);
	    hrecv_Zs.resize(Br.size3);
	    hrecv_Ze.resize(Br.size3);
	    omp_target_memcpy(hsend_Ze.data(), devptr_ro(Bs.Ze_data), sizeof(double)*Bs.size3, 0, 0, host, dev);
	    omp_target_memcpy(hsend_Zs.data(), devptr_ro(Bs.Zs_data), sizeof(double)*Bs.size3, 0, 0, host, dev);

	    double* h_Bs_Ze = hsend_Ze.data();
	    double* h_Bs_Zs = hsend_Zs.data();
	    double* h_Br_Zs = hrecv_Zs.data();
	    double* h_Br_Ze = hrecv_Ze.data();

    if (n3m != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Zs, Br.size3, MPI_DOUBLE, n3m, 3100, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Ze, Bs.size3, MPI_DOUBLE, n3m, 3200, comm3d, &req[nreq++]);
    } else {
      if (boundary_zin == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=js; j<=je; j++)
              for (int i=is; i<=ie; i++)
                Br.Zs(n,k,j,i) = Bs.Ze(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Zs(nve3,k,j,i) = -Br.Zs(nve3,k,j,i);
      } else if (boundary_zin == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=js; j<=je; j++)
              for (int i=is; i<=ie; i++)
                Br.Zs(n,k,j,i) = Bs.Ze(n,0,j,i);
      }
    }

    if (n3p != MPI_PROC_NULL) {
      rc = MPI_Irecv(h_Br_Ze, Br.size3, MPI_DOUBLE, n3p, 3200, comm3d, &req[nreq++]);
      rc = MPI_Isend(h_Bs_Zs, Bs.size3, MPI_DOUBLE, n3p, 3100, comm3d, &req[nreq++]);
    } else {
      if (boundary_zout == reflection) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=js; j<=je; j++)
              for (int i=is; i<=ie; i++)
                Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1-k,j,i);
#pragma omp target teams distribute parallel for collapse(3)
        for (int k=0; k<ngh; k++)
          for (int j=js; j<=je; j++)
            for (int i=is; i<=ie; i++)
              Br.Ze(nve3,k,j,i) = -Br.Ze(nve3,k,j,i);
      } else if (boundary_zout == outflow) {
#pragma omp target teams distribute parallel for collapse(4)
        for (int n=0; n<nprim; n++)
          for (int k=0; k<ngh; k++)
            for (int j=js; j<=je; j++)
              for (int i=is; i<=ie; i++)
                Br.Ze(n,k,j,i) = Bs.Zs(n,ngh-1,j,i);
      }
    }
  }

	  if (nreq != 0) {
	    MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);
	    // Copy host-received buffers -> device buffers (only for directions that used MPI)
	    if (ntiles[dir1] != 1) {
	      omp_target_memcpy(devptr_rw(Br.Xs_data), hrecv_Xs.data(), sizeof(double)*Br.size1, 0, 0, dev, host);
	      omp_target_memcpy(devptr_rw(Br.Xe_data), hrecv_Xe.data(), sizeof(double)*Br.size1, 0, 0, dev, host);
	    }
	    if (ntiles[dir2] != 1) {
	      omp_target_memcpy(devptr_rw(Br.Ys_data), hrecv_Ys.data(), sizeof(double)*Br.size2, 0, 0, dev, host);
	      omp_target_memcpy(devptr_rw(Br.Ye_data), hrecv_Ye.data(), sizeof(double)*Br.size2, 0, 0, dev, host);
	    }
	    if (ntiles[dir3] != 1) {
	      omp_target_memcpy(devptr_rw(Br.Zs_data), hrecv_Zs.data(), sizeof(double)*Br.size3, 0, 0, dev, host);
	      omp_target_memcpy(devptr_rw(Br.Ze_data), hrecv_Ze.data(), sizeof(double)*Br.size3, 0, 0, dev, host);
	    }
	  }
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
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  Bs.Xs(n,k,j,i) = P(n,k,j,ie-ngh+1+i);
	  Bs.Xe(n,k,j,i) = P(n,k,j,is      +i);
	  //Bs.Xs_data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+ie-ngh+1+i];
	  //Bs.Xe_data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+is      +i];
  }
  
  // y-direction
  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  Bs.Ys(n,k,j,i) = P(n,k,je-ngh+1+j,i);
	  Bs.Ye(n,k,j,i) = P(n,k,js      +j,i);
	  //Bs.Ys_data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+je-ngh+1+j)*itot+i];
	  //Bs.Ye_data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+js      +j)*itot+i];
  }
  
  // z-direction
  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  Bs.Zs(n,k,j,i) = P(n,ke-ngh+1+k,j,i);
	  Bs.Ze(n,k,j,i) = P(n,ks      +k,j,i);
	  //Bs.Zs_data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ke-ngh+1+k)*jtot+j)*itot+i];
	  //Bs.Ze_data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ks      +k)*jtot+j)*itot+i];
  }

  SendRecvBoundary(Bs,Br);

  
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xe|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  P(n,k,j,is-ngh+i) = Br.Xs(n,k,j,i);
	  P(n,k,j,ie+1  +i) = Br.Xe(n,k,j,i);
	  //P.data[((n*ktot+k)*jtot+j)*itot+is-ngh+i] = Xs.data[((n*ktot+k)*jtot+j)*ngh+i];
	  //P.data[((n*ktot+k)*jtot+j)*itot+ie+1  +i] = Xe.data[((n*ktot+k)*jtot+j)*ngh+i];
  }

  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,k,js-ngh+j,i) = Br.Ys(n,k,j,i);
	  P(n,k,je+1  +j,i) = Br.Ye(n,k,j,i);
	  //P.data[((n*ktot+k)*jtot+js-ngh+j)*itot+i] = Ys.data[((n*ktot+k)*ngh+j)*itot+i];
	  //P.data[((n*ktot+k)*jtot+je+1  +j)*itot+i] = Ye.data[((n*ktot+k)*ngh+j)*itot+i];
  }

  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,ks-ngh+k,j,i) = Br.Zs(n,k,j,i);
	  P(n,ke+1  +k,j,i) = Br.Ze(n,k,j,i);
	  //P.data[((n*ktot+ks-ngh+k)*jtot+j)*itot+i] = Zs.data[((n*ngh+k)*jtot+j)*itot+i];
	  //P.data[((n*ktot+ke+1  +k)*jtot+j)*itot+i] = Ze.data[((n*ngh+k)*jtot+j)*itot+i];
  }
};

};// end namespace
