/**
 * @file main.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include <algorithm>
#include <chrono>

#include "resolution.hpp"
#include "hydro.hpp"
#include "boundary.hpp"

#include "main.hpp"
#include "mpi_config.hpp"
#include "output.hpp"


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

#pragma omp target update to ( G.x1a_data[0:G.n1],G.x1b_data[0:G.n1])
#pragma omp target update to ( G.x2a_data[0:G.n2],G.x2b_data[0:G.n2])
#pragma omp target update to ( G.x3a_data[0:G.n3],G.x3b_data[0:G.n3])

  
}
static void GenerateProblem(hydflux_mod::GridArray<double>& G,hydflux_mod::FieldArray<double>& P,hydflux_mod::FieldArray<double>& U) {
  using namespace resolution_mod;
  using namespace hydflux_mod;

  double eint = 1.0e0;
  double denc = 1.0e0;
  csiso = sqrt(eint/denc);
  chg = 0.0e0;
#pragma omp target update to ( csiso,chg)

  double pres = denc*csiso*csiso;

  double xc = 0.5*(x1max + x1min);
  double yc = 0.5*(x2max + x2min);
  double zc = 0.5*(x3max + x3min);
  double sigma2x = pow((0.1e0*(x1max-x1min)),2);
  double sigma2y = pow((0.1e0*(x2max-x2min)),2);
  double sigma2z = pow((0.1e0*(x3max-x3min)),2);
  
  //printf("cs=%e",csiso);
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double x = G.x1b(i);
	double y = G.x2b(j);
	double z = G.x3b(k);
	  
	P(nden,k,j,i) = denc;
	P(nve1,k,j,i) = 0.3e0;
	P(nve2,k,j,i) = 0.3e0;
	P(nve3,k,j,i) = 0.3e0;
	P(nene,k,j,i)  = 1.0e6*exp(-( pow(x-xc,2)/sigma2x
				     +pow(y-yc,2)/sigma2y
	                             +pow(z-zc,2)/sigma2z )); //specific internel energy	
	P(npre,k,j,i)  =pres;
	P(ncsp,k,j,i) = csiso;
	
	P(nbm1,k,j,i) = 0.0e0;
	P(nbm2,k,j,i) = 0.0e0;
	P(nbm3,k,j,i) = 0.0e0;
	P(nbps,k,j,i) = 0.0e0;
    };
  
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    U(mden,k,j,i) = P(nden,k,j,i);
    U(mrv1,k,j,i) = P(nden,k,j,i)*P(nve1,k,j,i);
    U(mrv2,k,j,i) = P(nden,k,j,i)*P(nve2,k,j,i);
    U(mrv3,k,j,i) = P(nden,k,j,i)*P(nve3,k,j,i);
    double ekin = 0.5*P(nden,k,j,i)*( P(nve1,k,j,i)*P(nve1,k,j,i)
	        		     +P(nve2,k,j,i)*P(nve2,k,j,i)
				     +P(nve3,k,j,i)*P(nve3,k,j,i));
    double emag = 0.5              *( P(nbm1,k,j,i)*P(nbm1,k,j,i)
	        		     +P(nbm2,k,j,i)*P(nbm2,k,j,i)
				     +P(nbm3,k,j,i)*P(nbm3,k,j,i));
     U(meto,k,j,i) = P(nene,k,j,i)*P(nden,k,j,i) + ekin + emag;
     
     U(mbm1,k,j,i) = P(nbm1,k,j,i);
     U(mbm2,k,j,i) = P(nbm2,k,j,i);
     U(mbm3,k,j,i) = P(nbm3,k,j,i);
     U(mbps,k,j,i) = P(nbps,k,j,i);
  };
#pragma omp target update to ( U.data[0: U.size])
#pragma omp target update to ( P.data[0: P.size])

}



int main() {
  
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace boundary_mod;
  using namespace mpi_config_mod;
  // Follow Fortran main.f90 / output.f90 style flags
  const bool nooutput = false;
  const bool forceoutput = true;   // force output at start/end
  const bool usualoutput = false;  // regular output (subject to dtout)

	// Like Fortran (output.f90): logical,parameter :: binaryout=.true.
	// true  -> MPI binary output
	// false -> ASCII output (unf%05d.dat)
	const bool binaryout = false;

  periodic[dir1] = 1;
  periodic[dir2] = 1;
  periodic[dir3] = 1;
  ntiles[dir1]   = 1;
  ntiles[dir2]   = 2;
  ntiles[dir3]   = 1;
  InitializeMPI();
	// Configure I/O mode
  mpi_dataio_mod::binaryout = binaryout;
  
  if(myid_w == 0) printf("setup grids and fields\n");
  
  AllocateHydroVariables(G,U,Fx,Fy,Fz,P);
 
  AllocateBoundaryVariables(Bs,Br);
  
  if (myid_w == 0) printf("grid size for x y z = %i %i %i\n",ngrid1*ntiles[dir1],ngrid2*ntiles[dir2],ngrid3*ntiles[dir3]);
  
  GenerateGrid(G);
  GenerateProblem(G,P,U);
  // Force output at the initial state (Fortran: call Output(forceoutput))
  Output(forceoutput);

  if (myid_w == 0) printf("entering main loop\n");
  int step = 0;
  auto time_beg = std::chrono::high_resolution_clock::now();

  for (step=0;step<stepmax;step++){
    ControlTimestep(G); 
    if (myid_w==0 && step%300 ==0 && !nooutput) printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    //printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    SetBoundaryCondition(P,Bs,Br);
    EvaluateCh();
    GetNumericalFlux1(G,P,Fx);
    GetNumericalFlux2(G,P,Fy);
    GetNumericalFlux3(G,P,Fz);
    UpdateConservU(G,Fx,Fy,Fz,U);
    DampPsi(G,U);
    UpdatePrimitvP(U,P);

    time_sim += dt;
    //printf("dt=%e\n",dt);
    if (!nooutput) Output(usualoutput);
    //if (!nooutput) Output1D(usualoutput);

    if(time_sim > time_max) break;
    
  }

  //DeallocateHydroVariables(U,Fx,Fy,Fz,P);
  //DeallocateBoundaryVariables(Xs,Xe,Ys,Ye,Zs,Ze);
  
  auto time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = time_end - time_beg;
  if (myid_w == 0) printf("exiting main loop time=%e, step=%i\n",time_sim,step);
  if (myid_w == 0) printf("sim time [s]: %e\n", elapsed.count());
  if (myid_w == 0) printf("time/count/cell : %e\n", elapsed.count()/(ngrid1*ngrid2*ngrid3)/stepmax);

  // Force final output (Fortran: is_final=.true.; call Output(forceoutput))
  Output(forceoutput);
  //if (!nooutput) Output1D(forceoutput);
  
  printf("program has been finished\n");
  return 0;
}
