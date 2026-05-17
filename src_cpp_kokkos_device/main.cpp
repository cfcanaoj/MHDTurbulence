/**
 * @file main.cpp
 * @brief 
 * @author Keigo Nitadori
 * @date 2026-05-17
*/
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include <algorithm>
#include <chrono>
#include <random>

#include <Kokkos_Core.hpp>
#include <omp.h>

#include "config.hpp"
#include "mpi_config.hpp"
#include "mhd.hpp"
#include "boundary.hpp"
#include "output.hpp"
#include "main.hpp"

static void RealTimeAnalysis(const GridArray<double>& G, const FieldArray<double>& P) {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;

  double mix = 0.0;
  double avevy = 0.0;
  double vol = 0.0;

#if 0
// #pragma omp target teams distribute parallel for collapse(3) reduction(+:vol,mix,avevy)
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        const double dv = (G.x1a(i+1)-G.x1a(i)) * (G.x2a(j+1)-G.x2a(j)) * (G.x3a(k+1)-G.x3a(k));
        vol += dv;
        mix += P(nst,k,j,i) * (1.0 - P(nst,k,j,i)) * dv;
        avevy += P(nve2,k,j,i) * P(nve2,k,j,i) * dv;
      }
    }
  }
#else
	Kokkos::parallel_reduce("RealTimeAnalysis",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+1, ke+1}),
	KOKKOS_LAMBDA(const int i, const int j, const int k, double &lmix, double &lavevy, double &lvol) {
		const double dv = (G.dev_x1a(i+1)-G.dev_x1a(i)) * (G.dev_x2a(j+1)-G.dev_x2a(j)) * (G.dev_x3a(k+1)-G.dev_x3a(k));
		lvol += dv;
		lmix += P.dev(nst,k,j,i) * (1.0 - P.dev(nst,k,j,i)) * dv;
		lavevy += P.dev(nve2,k,j,i) * P.dev(nve2,k,j,i) * dv;
	},
	mix, avevy, vol);
#endif

  double local[3] = {vol, mix, avevy};
  double global[3] = {0.0, 0.0, 0.0};
  GetMPIsum(3, local, global);
  vol = global[0];
  mix = global[1] / vol;
  avevy = std::sqrt(global[2] / vol);

  if (myid_w == 0) {
    static bool is_inited = false;
    static FILE* fp = nullptr;
    constexpr double Amp = 1.2e-3;
    constexpr double Gamma = 1.49;
    if (!is_inited) {
      fp = std::fopen("t-prof.csv", "w");
      if (fp != nullptr) {
        std::fprintf(fp, "# Gamma= %24.16E\n", Gamma);
        std::fprintf(fp, "# 1:time 2:mix 3:v_y 4:exp(2Gamma*t)\n");
      }
      is_inited = true;
    }
    if (fp != nullptr) {
      std::fprintf(fp, " %24.16E %24.16E %24.16E %24.16E\n", time_sim, mix, avevy, Amp*std::exp(Gamma*time_sim));
      std::fflush(fp);
    }
  }
}

static void GenerateGrid(hydflux_mod::GridArray<double>& G) {
  using namespace resolution_mod;
  using namespace mpi_config_mod;
  
  double x1minloc,x1maxloc;
  double x2minloc,x2maxloc;
  double x3minloc,x3maxloc;
  double dx1,dx2,dx3;
   
  x1minloc = x1min + (x1max-x1min)/ntiles[dir1]* coords[dir1];
  x1maxloc = x1min + (x1max-x1min)/ntiles[dir1]*(coords[dir1]+1);
  
  dx1 = (x1maxloc-x1minloc)/double(ngrid1);
  for(int i=is-ngh;i<= ie+ngh+1;i++){
    G.x1a(i) = dx1*(i-(ngh+1))+x1minloc;
  }
  for(int i=is-ngh;i<= ie+ngh;i++){
    G.x1b(i) = 0.5e0*(G.x1a(i+1)+G.x1a(i));
  }

  x2minloc = x2min + (x2max-x2min)/ntiles[dir2]* coords[dir2];
  x2maxloc = x2min + (x2max-x2min)/ntiles[dir2]*(coords[dir2]+1);
  dx2=(x2maxloc-x2minloc)/double(ngrid2);
  for(int j=js-ngh;j<= je+ngh+1;j++){
    G.x2a(j) = dx2*(j-(ngh+1))+x2minloc;
  }
  for(int j=js-ngh;j<= je+ngh;j++){
    G.x2b(j) = 0.5e0*(G.x2a(j+1)+G.x2a(j));
  }

  x3minloc = x3min + (x3max-x3min)/ntiles[dir3]* coords[dir3];
  x3maxloc = x3min + (x3max-x3min)/ntiles[dir3]*(coords[dir3]+1);
  dx3=(x3maxloc-x3minloc)/double(ngrid3);
  for(int k=ks-ngh;k<= ke+ngh+1;k++){
    G.x3a(k) = dx3*(k-(ngh+1))+x3minloc;
  }
  for(int k=ks-ngh;k<= ke+ngh;k++){
    G.x3b(k) = 0.5e0*(G.x3a(k+1)+G.x3a(k));
  }
}
static void GenerateProblem(hydflux_mod::GridArray<double>& G,hydflux_mod::FieldArray<double>& P,hydflux_mod::FieldArray<double>& U) {
  // Kelvin-Helmholtz initial condition
  using namespace resolution_mod;
  using namespace mpi_config_mod;
  using namespace hydflux_mod;

  const double pi  = std::acos(-1.0);

  // Parameters (same as Fortran)
  const double rho1 = 1.0;
  const double rho2 = 1.0;
  const double dv   = 2.0;
  const double wid  = 0.05;
  const double sig  = 0.2;
  const double B0   = std::sqrt(2.0/3.0);
  (void)rho1; (void)rho2; (void)B0; // currently symmetric / B0 unused in reference file

  // Random perturbation strength (Fortran rrv)
  const double rrv = 1.0e-2;

  // Isothermal: choose csiso so that p = rho*csiso^2 matches Fortran p=1 when rho=1
  //csiso = 1.0;
  //#pragma omp target update to (csiso)
  
  chg   = 0.0;

  // Base profile (fill including ghost cells to avoid uninitialized values)
  for (int k = ks-ngh; k <= ke+ngh; ++k)
    for (int j = js-ngh; j <= je+ngh; ++j)
      for (int i = is-ngh; i <= ie+ngh; ++i) {
        const double x = G.x1b(i);
        const double y = G.x2b(j);

        const double tanh_p = std::tanh((y + 0.5)/wid);
        const double tanh_m = std::tanh((y - 0.5)/wid);

        P(nden,k,j,i) = 1.0;
        P(npre,k,j,i) = 1.0;

        P(nve1,k,j,i) = 0.5 * dv * (tanh_p - tanh_m - 1.0);
        P(nve2,k,j,i) = 1.0e-3 * std::sin(2.0*pi*x) *
                        ( std::exp(- (y + 0.5)*(y + 0.5)/(sig*sig))
                        + std::exp(- (y - 0.5)*(y - 0.5)/(sig*sig)) );
        P(nve3,k,j,i) = 0.0;

        P(nbm1,k,j,i) = 0.0;
        P(nbm2,k,j,i) = 0.0;
        P(nbm3,k,j,i) = 0.0;
        P(nbps,k,j,i) = 0.0;

        // specific internal energy is not used in the isothermal closure, but keep it consistent
        P(nene,k,j,i) = P(npre,k,j,i) /(gam-1.0);
	P(ncsp,k,j,i) = sqrt(gam*P(npre,k,j,i)/P(nden,k,j,i));
	  
	  //P(ncsp,k,j,i) = csiso;
        // Composition (Fortran: Xcomp(1)=0.5*(tanh((y+0.5)/wid)-tanh((y-0.5)/wid)))
        for (int n=0; n<ncomp; ++n) {
          (void)n; // for future multi-comp extension
          P(nst, k,j,i) = 0.5 * (tanh_p - tanh_m);
        }
      }

  // Random perturbation on v1 (Fortran: one random number per (j,k), applied to all i)
  // Fortran uses random_seed with seed(1)=1, seed(2)=1 + myid_w*in*jn*kn.
  // Here we emulate with a deterministic per-rank seed.
  std::mt19937_64 rng(static_cast<unsigned long long>(1 + myid_w * 1000003ULL));
  std::uniform_real_distribution<double> uni(0.0, 1.0);

  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++) {
      //      const double r = uni(rng);
      // const double dv_rand = dv * rrv * (r - 0.5);
      for (int i = is; i <= ie; ++i) {
        //P(nve1,k,j,i) += dv_rand;
      }
    }

  namespace res = resolution_mod; 
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++) {
      const double y = G.x2b(j);
      const double z = G.x3b(k);
      const double dv_harm = dv * rrv * std::cos(6*2*pi*(y - res::x2min)/(res::x2max - res::x2min))
                                      * std::cos(6*2*pi*(z - res::x3min)/(res::x3max - res::x3min));
      for (int i = is; i <= ie; i++) {
        P(nve1,k,j,i) += dv_harm;
      }
    }

  
  // Build conserved variables U from primitives P
  for (int k = ks-ngh; k <= ke+ngh; ++k)
    for (int j = js-ngh; j <= je+ngh; ++j)
      for (int i = is-ngh; i <= ie+ngh; ++i) {
        const double rho = P(nden,k,j,i);
        U(mden,k,j,i) = rho;
        U(mrv1,k,j,i) = rho * P(nve1,k,j,i);
        U(mrv2,k,j,i) = rho * P(nve2,k,j,i);
        U(mrv3,k,j,i) = rho * P(nve3,k,j,i);

        const double ekin = 0.5 * rho * ( P(nve1,k,j,i)*P(nve1,k,j,i)
                                       + P(nve2,k,j,i)*P(nve2,k,j,i)
                                       + P(nve3,k,j,i)*P(nve3,k,j,i) );
        const double emag = 0.5 * ( P(nbm1,k,j,i)*P(nbm1,k,j,i)
                                  + P(nbm2,k,j,i)*P(nbm2,k,j,i)
                                  + P(nbm3,k,j,i)*P(nbm3,k,j,i) );

        U(meto,k,j,i) = rho * P(nene,k,j,i) + ekin + emag;

        U(mbm1,k,j,i) = P(nbm1,k,j,i);
        U(mbm2,k,j,i) = P(nbm2,k,j,i);
        U(mbm3,k,j,i) = P(nbm3,k,j,i);
        U(mbps,k,j,i) = P(nbps,k,j,i);

        for (int n=0; n<ncomp; ++n) {
          U(mst+n,k,j,i) = rho * P(nst+n,k,j,i);
        }
      }
}



int main(int argc, char **argv) {
  
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace boundary_mod;
  using namespace mpi_config_mod;
  const bool forceoutput = true;   // force output at start/end
  const bool usualoutput = false;  // regular output (subject to dtout)

  InitializeMPI();
  Kokkos::initialize(argc, argv);

#ifdef _OPENMP
  int nth = -1;
#pragma omp parallel
  nth = omp_get_num_threads();
  printf("%d: %d threads\n", myid_w, nth);
#endif

  
  if(myid_w == 0) printf("setup grids and fields\n");
  
  AllocateHydroVariables(G,U,Fx,Fy,Fz,P);
 
  mpi_config_mod::gpu_aware = true;
  AllocateBoundaryVariables(Bs,Br);
  
  if (myid_w == 0) printf("grid size for x y z = %i %i %i\n",ngrid1*ntiles[dir1],ngrid2*ntiles[dir2],ngrid3*ntiles[dir3]);
  
  GenerateGrid(G);
  GenerateProblem(G,P,U);
  G.h2d();
  P.h2d();
  U.h2d();
  RealTimeAnalysis(G, P); // currently on the host
  // Force output at the initial state (Fortran: call Output(forceoutput))
  Output(forceoutput);

  if (myid_w == 0) printf("entering main loop\n");
  int step = 0;
  auto time_beg = std::chrono::high_resolution_clock::now();

  FILE *fpdt = nullptr;
  if(!myid_w) fpdt = fopen("dt.log", "w");

#define SIMPLE_PROFILING
  auto wtime = []{
#ifdef SIMPLE_PROFILING
	  Kokkos::fence();
	  MPI_Barrier(MPI_COMM_WORLD);
	  return MPI_Wtime();
#else
	  return 0.0;
#endif
  };
  double wt_dt=0.0, wt_bound=0.0, wt_ch=0.0, wt_flux=0.0, wt_upC=0.0, wt_dump=0.0, wt_upP=0.0;

  for (step=0;step<stepmax;step++){

    double t0 = wtime();
    ControlTimestep(G, P); 
    double t1 = wtime();

    // if (myid_w==0 && step%300 ==0 && ! config::benchmarkmode) printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    //if (myid_w==0 && step%10 ==0) printf("step=%i/%i time=%e dt=%e, %f%%\n",step,stepmax,time_sim,dt, time_sim/time_max*100.0);
    //printf("step=%i time=%e dt=%e\n",step,time_sim,dt);

    double t2 = wtime();
    SetBoundaryCondition(P,Bs,Br);
    double t3 = wtime();

    EvaluateCh(P);
    double t4 = wtime();
    /*
    if(fpdt && step<100){
	    fprintf(fpdt, "step=%d, time=%e, dt=%e, chg=%e\n", step, time_sim, dt, chg);
    }
    if(fpdt && step>=100){
	    fclose(fpdt);
	    fpdt = nullptr;
    }
    */

    double t5 = wtime();
    GetNumericalFluxD(1, G,P,Fx);
    GetNumericalFluxD(2, G,P,Fy);
    GetNumericalFluxD(3, G,P,Fz);
    double t6 = wtime();

    UpdateConservU(G,Fx,Fy,Fz,U);
    double t7 = wtime();

    DampPsi(G,U);
    double t8 = wtime();

    UpdatePrimitvP(U,P);
    double t9 = wtime();

    time_sim += dt;

    wt_dt    += t1-t0;
    wt_bound += t3-t2;
    wt_ch    += t4-t3;
    wt_flux  += t6-t5;
    wt_upC   += t7-t6;
    wt_dump  += t8-t7;
    wt_upP   += t9-t8;
    // printf("dt=%e\n",dt);
    if (! config::benchmarkmode) RealTimeAnalysis(G,P);
    if (! config::benchmarkmode) Output(usualoutput);
    //if (!nooutput) Output1D(usualoutput);

    if(time_sim > time_max) break;
    //if(time_sim > time_max || step >= 1000) break;
    
  }
  P.d2h();

  auto time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = time_end - time_beg;
  if (myid_w == 0) printf("exiting main loop time=%e, step=%i\n",time_sim,step);
  if (myid_w == 0) printf("sim time [s]: %e\n", elapsed.count());
  if (myid_w == 0) printf("time/count/cell : %e\n", elapsed.count()/(ngrid1*ngrid2*ngrid3*ntiles[dir1]*ntiles[dir2]*ntiles[dir3])/(step+1));
#ifdef SIMPLE_PROFILING
  if (myid_w == 0) printf("%e %e %e %e %e %e %e\n", wt_dt, wt_bound, wt_ch, wt_flux, wt_upC, wt_dump, wt_upP);
#endif

  // Force final output (Fortran: is_final=.true.; call Output(forceoutput))
  Output(forceoutput);
  RealTimeAnalysis(G, P);
  //if (!nooutput) Output1D(forceoutput);
  
  if (myid_w == 0) printf("program has been finished\n");

  DeallocateHydroVariables(G,U,Fx,Fy,Fz,P);
  DeallocateBoundaryVariables(Bs, Br);

  Kokkos::finalize();
  FinalizeMPI();

  return 0;
}
