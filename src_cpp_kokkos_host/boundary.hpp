/**
 * @file boundary.hpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include "config.hpp"
#include "mhd.hpp"
using namespace hydflux_mod;

namespace boundary_mod {

  inline constexpr int periodicb   = config::periodicb;
  inline constexpr int reflection  = config::reflection;
  inline constexpr int outflow     = config::outflow;

  // Per-face boundary selection (1:periodic, 2:reflection, 3:outflow)
  extern int boundary_xin, boundary_xout;
  extern int boundary_yin, boundary_yout;
  extern int boundary_zin, boundary_zout;

#define VIEW_BOUNDARY
  template <typename T>
  class BoundaryArray {
  public:
    T* Xs_data = nullptr;
    T* Xe_data = nullptr;
    T* Ys_data = nullptr;
    T* Ye_data = nullptr;
    T* Zs_data = nullptr;
    T* Ze_data = nullptr;
    int nv=0, ng=0, n3 = 0, n2 = 0, n1 = 0;
    size_t size1 = 0;
    size_t size2 = 0;
    size_t size3 = 0;
    
    BoundaryArray() = default;
    BoundaryArray(int nv_,int ng_,int n3_, int n2_, int n1_) {
      allocate(nv_,ng_,n3_,n2_,n1_);
    }
#ifndef VIEW_BOUNDARY
    void allocate(int nv_,int ng_,int n3_, int n2_, int n1_) {
      nv = nv_;ng = ng_; n3 = n3_; n2 = n2_; n1 = n1_;
      size1 = static_cast<size_t>(nv) * n3 * n2 * ng;
      Xs_data = static_cast<T*>(malloc(sizeof(T) * size1));
      Xe_data = static_cast<T*>(malloc(sizeof(T) * size1));
      size2 = static_cast<size_t>(nv) * n3 * ng * n1;
      Ys_data = static_cast<T*>(malloc(sizeof(T) * size2));
      Ye_data = static_cast<T*>(malloc(sizeof(T) * size2));
      size3 = static_cast<size_t>(nv) * ng * n2 * n1;
      Zs_data = static_cast<T*>(malloc(sizeof(T) * size3));
      Ze_data = static_cast<T*>(malloc(sizeof(T) * size3));
    }
    
    inline       T& Xs(int n, int k, int j, int i)       noexcept { return Xs_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline const T& Xs(int n, int k, int j, int i) const noexcept { return Xs_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline       T& Xe(int n, int k, int j, int i)       noexcept { return Xe_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline const T& Xe(int n, int k, int j, int i) const noexcept { return Xe_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline       T& Ys(int n, int k, int j, int i)       noexcept { return Ys_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline const T& Ys(int n, int k, int j, int i) const noexcept { return Ys_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline       T& Ye(int n, int k, int j, int i)       noexcept { return Ye_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline const T& Ye(int n, int k, int j, int i) const noexcept { return Ye_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline       T& Zs(int n, int k, int j, int i)       noexcept { return Zs_data[((n*ng + k)*n2 + j)*n1 + i]; }
    inline const T& Zs(int n, int k, int j, int i) const noexcept { return Zs_data[((n*ng + k)*n2 + j)*n1 + i]; }
    inline       T& Ze(int n, int k, int j, int i)       noexcept { return Ze_data[((n*ng + k)*n2 + j)*n1 + i]; }
    inline const T& Ze(int n, int k, int j, int i) const noexcept { return Ze_data[((n*ng + k)*n2 + j)*n1 + i]; }
#else
    using DView = Kokkos::View<T****, Kokkos::LayoutLeft>;
    using HView = typename DView::HostMirror;
    DView d_Xs, d_Xe, d_Ys, d_Ye, d_Zs, d_Ze;
    HView h_Xs, h_Xe, h_Ys, h_Ye, h_Zs, h_Ze;

    void allocate(int nv_,int ng_,int n3_, int n2_, int n1_) {
      nv = nv_;ng = ng_; n3 = n3_; n2 = n2_; n1 = n1_;

      d_Xs = DView("Xs", ng, nv, n2, n3);
      d_Xe = DView("Xe", ng, nv, n2, n3);
      d_Ys = DView("Ys", n1, nv, ng, n3);
      d_Ye = DView("Ye", n1, nv, ng, n3);
      d_Zs = DView("Zs", n1, nv, n2, ng);
      d_Ze = DView("Ze", n1, nv, n2, ng);

      h_Xs = Kokkos::create_mirror_view(d_Xs);
      h_Xe = Kokkos::create_mirror_view(d_Xe);
      h_Ys = Kokkos::create_mirror_view(d_Ys);
      h_Ye = Kokkos::create_mirror_view(d_Ye);
      h_Zs = Kokkos::create_mirror_view(d_Zs);
      h_Ze = Kokkos::create_mirror_view(d_Ze);
      
      Xs_data = h_Xs.data();
      Xe_data = h_Xe.data();
      Ys_data = h_Ys.data();
      Ye_data = h_Ye.data();
      Zs_data = h_Zs.data();
      Ze_data = h_Ze.data();

      size1 = d_Xs.size();
      size2 = d_Ys.size();
      size3 = d_Zs.size();

#if 0
      printf("(%p:%p),(%p:%p),(%p:%p)\n",
	      Xs_data, Xe_data, Ys_data, Ye_data, Zs_data, Ze_data);
#endif
    }

    // copy captureされても書き込めるようにconstは外した
    inline       T& Xs(int n, int k, int j, int i)       noexcept { return h_Xs(i, n, j, k); }
    inline       T& Xs(int n, int k, int j, int i) const noexcept { return h_Xs(i, n, j, k); }
    inline       T& Xe(int n, int k, int j, int i)       noexcept { return h_Xe(i, n, j, k); }
    inline       T& Xe(int n, int k, int j, int i) const noexcept { return h_Xe(i, n, j, k); }
    inline       T& Ys(int n, int k, int j, int i)       noexcept { return h_Ys(i, n, j, k); }
    inline       T& Ys(int n, int k, int j, int i) const noexcept { return h_Ys(i, n, j, k); }
    inline       T& Ye(int n, int k, int j, int i)       noexcept { return h_Ye(i, n, j, k); }
    inline       T& Ye(int n, int k, int j, int i) const noexcept { return h_Ye(i, n, j, k); }
    inline       T& Zs(int n, int k, int j, int i)       noexcept { return h_Zs(i, n, j, k); }
    inline       T& Zs(int n, int k, int j, int i) const noexcept { return h_Zs(i, n, j, k); }
    inline       T& Ze(int n, int k, int j, int i)       noexcept { return h_Ze(i, n, j, k); }
    inline       T& Ze(int n, int k, int j, int i) const noexcept { return h_Ze(i, n, j, k); }
#endif

    void deallocate(){
	    d_Xs = d_Xe = d_Ys = d_Ye = d_Zs = d_Ze = DView();
	    h_Xs = h_Xe = h_Ys = h_Ye = h_Zs = h_Ze = HView();
    }

    void d2h_Xs(){ Kokkos::deep_copy(h_Xs, d_Xs); }
    void d2h_Xe(){ Kokkos::deep_copy(h_Xe, d_Xe); }

    void h2d_Xs(){ Kokkos::deep_copy(d_Xs, h_Xs); }
    void h2d_Xe(){ Kokkos::deep_copy(d_Xe, h_Xe); }
  };

  extern BoundaryArray<double> Bs,Br; 

  void AllocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
  void DeallocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
  
  void SetBoundaryCondition(FieldArray<double>& P,BoundaryArray<double>& Bs,BoundaryArray<double>& Br);

  void SendRecvBoundary(const BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
};

#endif
