/**
 * @file mhd.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#ifndef MHD_HPP_
#define MHD_HPP_

#include <cmath>
#include <cstdio>
#include <algorithm>

#include <Kokkos_Core.hpp>

namespace hydflux_mod {
  using index_t = int;     // 必要なら 32/64 を切り替え
  using size_t  = std::size_t;
  
#define GRID_VIEW
  template <typename T>
  class GridArray {
  public:
    T* x1a_data = nullptr;
    T* x1b_data = nullptr;
    T* x2a_data = nullptr;
    T* x2b_data = nullptr;
    T* x3a_data = nullptr;
    T* x3b_data = nullptr;
    int n3 = 0, n2 = 0, n1 = 0;
    size_t size = 0;
    GridArray() = default;
    GridArray(int n3_, int n2_, int n1_) {
      allocate(n3_,n2_,n1_);
    }
#ifndef GRID_VIEW
    void allocate(int n3_, int n2_, int n1_) {
      n3 = n3_; n2 = n2_; n1 = n1_;
      size = static_cast<size_t>(n1);
      x1a_data = static_cast<T*>(malloc(sizeof(T) * size));
      x1b_data = static_cast<T*>(malloc(sizeof(T) * size));
      size = static_cast<size_t>(n2);
      x2a_data = static_cast<T*>(malloc(sizeof(T) * size));
      x2b_data = static_cast<T*>(malloc(sizeof(T) * size));
      size = static_cast<size_t>(n3);
      x3a_data = static_cast<T*>(malloc(sizeof(T) * size));
      x3b_data = static_cast<T*>(malloc(sizeof(T) * size));
    }
    
    inline T& x1a(int i)       noexcept { return x1a_data[i]; }
    inline T& x1b(int i)       noexcept { return x1b_data[i]; }
    // const
    inline const T& x1a(int i) const noexcept { return x1a_data[i]; }
    inline const T& x1b(int i) const noexcept { return x1b_data[i]; }

    inline T& x2a(int j)       noexcept { return x2a_data[j]; }
    inline T& x2b(int j)       noexcept { return x2b_data[j]; }
    // const
    inline const T& x2a(int j) const noexcept { return x2a_data[j]; }
    inline const T& x2b(int j) const noexcept { return x2b_data[j]; }
  
    inline T& x3a(int k)       noexcept { return x3a_data[k]; }
    inline T& x3b(int k)       noexcept { return x3b_data[k]; }
    // const
    inline const T& x3a(int k) const noexcept { return x3a_data[k]; }
    inline const T& x3b(int k) const noexcept { return x3b_data[k]; }
#else
    using DView = Kokkos::View<T**, Kokkos::LayoutLeft>;
    using HView = typename DView::HostMirror;

    DView d1, d2, d3;
    HView h1, h2, h3;

    void allocate(int n3_, int n2_, int n1_) {
      n3 = n3_; n2 = n2_; n1 = n1_;

      d1 = DView("Grid_x", n1, 2);
      d2 = DView("Grid_x", n2, 2);
      d3 = DView("Grid_x", n3, 2);

      h1 = Kokkos::create_mirror_view(d1);
      h2 = Kokkos::create_mirror_view(d2);
      h3 = Kokkos::create_mirror_view(d3);
    }

    void deallocate() {
	    h1 = h2 = h3 = HView();
	    d1 = d2 = d3 = DView();
    }

    T& x1a(int i) noexcept { return h1(i, 0); }
    T& x1b(int i) noexcept { return h1(i, 1); }
    T& x2a(int j) noexcept { return h2(j, 0); }
    T& x2b(int j) noexcept { return h2(j, 1); }
    T& x3a(int k) noexcept { return h3(k, 0); }
    T& x3b(int k) noexcept { return h3(k, 1); }

    const T& x1a(int i) const noexcept { return h1(i, 0); }
    const T& x1b(int i) const noexcept { return h1(i, 1); }
    const T& x2a(int j) const noexcept { return h2(j, 0); }
    const T& x2b(int j) const noexcept { return h2(j, 1); }
    const T& x3a(int k) const noexcept { return h3(k, 0); }
    const T& x3b(int k) const noexcept { return h3(k, 1); }

    void h2d(){
	    Kokkos::deep_copy(d1, h1);
	    Kokkos::deep_copy(d2, h2);
	    Kokkos::deep_copy(d3, h3);
    }
#endif

  };

#define VIEW_LAYOUT_INJK 
// #define VIEW_LAYOUT_IJKN 
  template <typename T>
  class FieldArray {
  public:
#if (!defined VIEW_LAYOUT_INJK) && (!defined VIEW_LAYOUT_IJKN)
// #warning w/o KOKKOS VIEW
    T* data = nullptr;
    int nv = 0, n3 = 0, n2 = 0, n1 = 0;
    size_t size = 0;
    FieldArray() = default;
    FieldArray(int _nv, int _n3, int _n2, int _n1) {
      allocate(_nv,_n3,_n2,_n1);
    }
    
    void allocate(int _nv, int _n3, int _n2, int _n1) {
      nv = _nv; n3 = _n3; n2 = _n2; n1 = _n1;
      size = static_cast<size_t>(nv) * n3 * n2 * n1;
      data = static_cast<T*>(malloc(sizeof(T) * size));
  }
    
    void deallocate() {
      free(data);
      data = nullptr;
    }
    
    inline const T& operator()(int n, int k, int j, int i) const noexcept {
      return data[((n*n3 + k)*n2 + j)*n1 + i];
    }
    inline  T& operator()(int n, int k, int j, int i)  noexcept {
      return data[((n*n3 + k)*n2 + j)*n1 + i];
    }

#else
// #warning USE KOKKOS VIEW
    using DView = Kokkos::View<T****, Kokkos::LayoutLeft>;
    using HView = typename DView::HostMirror;

    DView d_view;
    HView h_view;
    T* data = nullptr;
    int nv = 0, n3 = 0, n2 = 0, n1 = 0;
    size_t size = 0;

    FieldArray() = default;

    FieldArray(int _nv, int _n3, int _n2, int _n1) {
	    allocate(_nv,_n3,_n2,_n1);
    }

    void allocate(int _nv, int _n3, int _n2, int _n1, const std::string &name="F") {
	    nv = _nv; n3 = _n3; n2 = _n2; n1 = _n1;
	    size = static_cast<size_t>(nv) * n3 * n2 * n1;
#ifdef VIEW_LAYOUT_INJK
	    d_view = DView(name, n1, nv, n2, n3);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    d_view = DView(name, n1, n2, n3, nv);
#endif
	    h_view = Kokkos::create_mirror_view(d_view);
	    data = h_view.data();
    }
    
    void deallocate() {
	    h_view = HView();
	    d_view = DView();
	    data = nullptr;
    }

    const T& operator()(int n, int k, int j, int i) const noexcept {
#ifdef VIEW_LAYOUT_INJK
	    return h_view(i, n, j, k);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    return h_view(i, j, k, n);
#endif
    }

    T& operator()(int n, int k, int j, int i)  noexcept {
#ifdef VIEW_LAYOUT_INJK
	    return h_view(i, n, j, k);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    return h_view(i, j, k, n);
#endif
    }

    const T& dev(int n, int k, int j, int i) const noexcept {
#ifdef VIEW_LAYOUT_INJK
	    return d_view(i, n, j, k);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    return d_view(i, j, k, n);
#endif
    }

    T& dev(int n, int k, int j, int i)  noexcept {
#ifdef VIEW_LAYOUT_INJK
	    return d_view(i, n, j, k);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    return d_view(i, j, k, n);
#endif
    }

    // copy captureでconstになったときでも書き込めるように
    T& href(int n, int k, int j, int i) const noexcept {
#ifdef VIEW_LAYOUT_INJK
	    return h_view(i, n, j, k);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    return h_view(i, j, k, n);
#endif
    }

    // 同上、デバイス版
    T& dref(int n, int k, int j, int i) const noexcept {
#ifdef VIEW_LAYOUT_INJK
	    return d_view(i, n, j, k);
#endif
#ifdef VIEW_LAYOUT_IJKN
	    return d_view(i, j, k, n);
#endif
    }


    void h2d(){
	    Kokkos::deep_copy(d_view, h_view);
    }
    void d2h(){
	    Kokkos::deep_copy(h_view, d_view);
    }
#endif // Kokkos view version?
  };

  inline constexpr int ncomp{1}; //! composition
  inline constexpr int nprim{11+ncomp}; //!
  inline constexpr int nden{0},nve1{1},nve2{2},nve3{3},nene{4},npre{5},ncsp{6},
                               nbm1{7},nbm2{8},nbm3{9},nbps{10},
                               nst{nprim-ncomp},ned{nprim-1};
  
  inline constexpr int mconsv{9+ncomp},madd{3}; // total number!
  inline constexpr int mudn{ 0},muvu{ 1},muvv{ 2},muvw{ 3},muet{ 4},// 0 <= index <= total-1 
                                mubu{ 5},mubv{ 6},mubw{ 7},mubp{ 8},
                       mfdn{mconsv},mfvu{mconsv+1},mfvv{mconsv+2},mfvw{mconsv+3},mfet{mconsv+4},
                                    mfbu{mconsv+5},mfbv{mconsv+6},mfbw{mconsv+7},mfbp{mconsv+8},
                                    must{  mconsv-ncomp},mued{  mconsv-1},
                                    mfst{2*mconsv-ncomp},mfed{2*mconsv-1},
                       mcsp{2*mconsv},mvel{2*mconsv+1},mpre{2*mconsv+2};
  
  inline constexpr int mden{ 0},mrv1{ 1},mrv2{ 2},mrv3{ 3},meto{ 4}, // 0 <= index <= total-1
                                mbm1{ 5},mbm2{ 6},mbm3{ 7},mbps{ 8},
                                mrvu{ 1},mrvv{ 2},mrvw{ 3},
                                mbmu{ 5},mbmv{ 6},mbmw{ 7},
                                mst{mconsv-ncomp},med{mconsv-1};
  extern GridArray<double> G;
  extern FieldArray<double> P; //! P(nprim ,ktot,jtot,itot)
  extern FieldArray<double> U; //! U(mconsv,ktot,jtot,itot)
  extern FieldArray<double> Fx,Fy,Fz;
  //extern double csiso;
  extern double gam;
  extern double chg;
  
  void AllocateHydroVariables(GridArray<double>& G,FieldArray<double>& U,FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P);
  void DeallocateHydroVariables(GridArray<double>& G,FieldArray<double>& U, FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P);
  void GetNumericalFlux1(const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double>& Fx);
  void GetNumericalFlux2(const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double>& Fy);
  void GetNumericalFlux3(const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double>& Fz);
  void GetNumericalFluxD(const int dir, const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double> &F);
  void UpdateConservU(const GridArray<double>& G,const FieldArray<double>& Fx,const FieldArray<double>& Fy,const FieldArray<double>& Fz,FieldArray<double>& U);
  void UpdatePrimitvP(const FieldArray<double>& U,FieldArray<double>& P);
  void ControlTimestep(const GridArray<double>& G);
  void EvaluateCh();
  void DampPsi(const GridArray<double>& G,FieldArray<double>& U);
  
};
#endif
