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
#include "mpi_dataio.hpp"


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
  // Kelvin–Helmholtz instability initial condition (ported from main.f90)
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;

  const double pi = std::acos(-1.0e0);
  const double rho1 = 1.0e0;
  const double rho2 = 1.0e0;
  const double dv   = 2.0e0;
  const double wid  = 0.05e0;
  const double sig  = 0.2e0;
  const double rrv  = 1.0e-2; // random perturbation strength (fraction of dv)

  // Isothermal: p = rho * cs^2. Choose cs=1 so p=1 for rho=1.
  const double eint = 1.0e0;
  const double d0   = 1.0e0;
  csiso = std::sqrt(eint/d0);
  chg   = 0.0e0;
#pragma omp target update to (csiso, chg)

  // deterministic hash-based pseudo-random in [0,1)
  auto u01 = [&](int a, int b, int c) -> double {
    unsigned long long x = 1469598103934665603ull;
    auto mix = [&](unsigned long long v){ x ^= v + 0x9e3779b97f4a7c15ull + (x<<6) + (x>>2); };
    mix((unsigned long long)a);
    mix((unsigned long long)b);
    mix((unsigned long long)c);
    // xorshift64*
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    unsigned long long y = x * 2685821657736338717ull;
    // take top 53 bits
    return (double)((y >> 11) & ((1ull<<53)-1)) / (double)(1ull<<53);
  };

  // fill primitives (including ghosts for periodic boundaries)
  for (int k=ks-ngh; k<=ke+ngh; ++k)
    for (int j=js-ngh; j<=je+ngh; ++j)
      for (int i=is-ngh; i<=ie+ngh; ++i) {
        const double x = G.x1b(i);
        const double y = G.x2b(j);

        const double rho = (y < 0.0e0) ? rho1 : rho2;
        const double pres = 1.0e0;

        // v1 profile (shear layers around y=±0.5)
        double v1 = 0.5e0 * dv * ( std::tanh((y + 0.5e0)/wid) - std::tanh((y - 0.5e0)/wid) - 1.0e0 );

        // v2 perturbation
        double v2 = 1.0e-3 * std::sin(2.0e0*pi*x) * ( std::exp(- (y + 0.5e0)*(y + 0.5e0)/(sig*sig))
                                                     + std::exp(- (y - 0.5e0)*(y - 0.5e0)/(sig*sig)) );
        double v3 = 0.0e0;

        // random perturbation to v1 (same random number for each (j,k) in Fortran; here make it depend on (j,k))
        const double r = u01(j, k, myid_w) - 0.5e0;
        v1 += dv * rrv * r;

        // passive scalar Xcomp: layer indicator (0..1)
        const double xc = 0.5e0*( std::tanh((y + 0.5e0)/wid) - std::tanh((y - 0.5e0)/wid) );

        P(nden,k,j,i) = rho;
        P(nve1,k,j,i) = v1;
        P(nve2,k,j,i) = v2;
        P(nve3,k,j,i) = v3;
        P(npre,k,j,i) = pres;
        P(ncsp,k,j,i) = csiso;
        P(nene,k,j,i) = eint/rho; // used only to build total energy

        P(nbm1,k,j,i) = 0.0e0;
        P(nbm2,k,j,i) = 0.0e0;
        P(nbm3,k,j,i) = 0.0e0;
        P(nbps,k,j,i) = 0.0e0;
        P(nxc ,k,j,i) = xc;
      }

  // build conserved
  for (int k=ks-ngh; k<=ke+ngh; ++k)
    for (int j=js-ngh; j<=je+ngh; ++j)
      for (int i=is-ngh; i<=ie+ngh; ++i) {
        U(mden,k,j,i) = P(nden,k,j,i);
        U(mrv1,k,j,i) = P(nden,k,j,i)*P(nve1,k,j,i);
        U(mrv2,k,j,i) = P(nden,k,j,i)*P(nve2,k,j,i);
        U(mrv3,k,j,i) = P(nden,k,j,i)*P(nve3,k,j,i);

        const double ekin = 0.5e0*P(nden,k,j,i)*( P(nve1,k,j,i)*P(nve1,k,j,i)
                                                 +P(nve2,k,j,i)*P(nve2,k,j,i)
                                                 +P(nve3,k,j,i)*P(nve3,k,j,i));
        const double emag = 0.5e0*( P(nbm1,k,j,i)*P(nbm1,k,j,i)
                                   +P(nbm2,k,j,i)*P(nbm2,k,j,i)
                                   +P(nbm3,k,j,i)*P(nbm3,k,j,i));
        U(meto,k,j,i) = P(nene,k,j,i)*P(nden,k,j,i) + ekin + emag;

        U(mbm1,k,j,i) = P(nbm1,k,j,i);
        U(mbm2,k,j,i) = P(nbm2,k,j,i);
        U(mbm3,k,j,i) = P(nbm3,k,j,i);
        U(mbps,k,j,i) = P(nbps,k,j,i);
        U(mxc ,k,j,i) = P(nden,k,j,i)*P(nxc ,k,j,i);
      }

#pragma omp target update to (U.data[0:U.size])
#pragma omp target update to (P.data[0:P.size])
}

static void GenerateProblem2(hydflux_mod::GridArray<double>& G,hydflux_mod::FieldArray<double>& P,hydflux_mod::FieldArray<double>& U) {
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


// ---- XMF writer for VisIt/ParaView (rank 0 only) ----
static void WriteXMF(int timeid, double time_sim) {
  using namespace mpi_config_mod;
  using namespace resolution_mod;
  namespace io = mpi_dataio_mod;
  if (myid_w != 0) return;

  char fname[256];
  std::snprintf(fname, sizeof(fname), "bindata/field%05d.xmf", timeid);
  FILE* fp = std::fopen(fname, "w");
  if (!fp) return;

  const int NX = io::ntotal[0];
  const int NY = io::ntotal[1];
  const int NZ = io::ntotal[2];

  const double dx = (x1max - x1min) / double(NX);
  const double dy = (x2max - x2min) / double(NY);
  const double dz = (x3max - x3min) / double(NZ);

  char datafile[256];
  std::snprintf(datafile, sizeof(datafile), "bindata/d3d%s.%05d", io::id, timeid);

  auto attr = [&](const char* name, int vid){
    std::fprintf(fp,
      "    <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
      "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n"
      "        <DataItem Dimensions=\"3 4\" Format=\"XML\">%d 0 0 0  1 1 1 1  1 %d %d %d</DataItem>\n"
      "        <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\">%s</DataItem>\n"
      "      </DataItem>\n"
      "    </Attribute>\n",
      name, NZ, NY, NX,
      vid, NZ, NY, NX,
      io::nvars, NZ, NY, NX,
      datafile);
  };

  std::fprintf(fp,
    "<?xml version=\"1.0\" ?>\n"
    "<Xdmf Version=\"3.0\">\n"
    "  <Domain>\n"
    "    <Grid Name=\"mesh\" GridType=\"Uniform\">\n"
    "      <Time Value=\"%.17g\"/>\n"
    "      <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n"
    "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
    "        <DataItem Dimensions=\"3\" Format=\"XML\">%.17g %.17g %.17g</DataItem>\n"
    "        <DataItem Dimensions=\"3\" Format=\"XML\">%.17g %.17g %.17g</DataItem>\n"
    "      </Geometry>\n",
    time_sim, NZ+1, NY+1, NX+1, x3min, x2min, x1min, dz, dy, dx);

  attr("d", 0);
  attr("v1", 1);
  attr("v2", 2);
  attr("v3", 3);
  attr("b1", 4);
  attr("b2", 5);
  attr("b3", 6);
  attr("psi", 7);
  attr("p", 8);
  attr("Xcomp", 9);

  std::fprintf(fp,
    "    </Grid>\n"
    "  </Domain>\n"
    "</Xdmf>\n");
  std::fclose(fp);
}

void Output1D(bool& forcedamp){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  static int index = 0;  
  static bool is_inited = false;
  const int dir1=0, dir2=1, dir3=2;
  if(!forcedamp && time_sim < time_out + dtout) return;

  if(myid_w==0) printf("output index=%i, time=%e \n",index,time_sim);

#pragma omp target update from (P.data[0:P.size])
  int ic,jc,kc;

  // assuming division is up to 2
  if(ntiles[dir1] == 1) {
    ic = int((is+ie)/2);
  }else{
    if(coords[dir1] ==0) ic = ie;
    if(coords[dir1] ==1) ic = is;
  }

  if(ntiles[dir2] == 1) {
    jc = int((js+je)/2);
  }else{
    if(coords[dir2] ==0) jc = je;
    if(coords[dir2] ==1) jc = js;
  }

  
  if(ntiles[dir3] == 1) {
    kc = int((ks+ke)/2);
  }else{
    if(coords[dir3] ==0) kc = ke;
    if(coords[dir3] ==1) kc = ks;
  }
  
  if (! is_inited){
    (void)system("mkdir -p snap");
    is_inited = true;
  }
  // for x
  char fname[256];
  std::snprintf(fname, sizeof(fname), "snap/snx_%02d_%05d.dat",myid_w, index);
  FILE* fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int i=is-2;i<=ie+2;i++){
    std::fprintf(fp, "%e %e %e \n", G.x1b(i),P(nden,kc,jc,i),P(nene,kc,jc,i) );
  }
  std::fclose(fp);
  
  // for y
  std::snprintf(fname, sizeof(fname), "snap/sny_%02d_%05d.dat",myid_w, index);
  fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int j=js-2;j<=je+2;j++){
    std::fprintf(fp, "%e %e %e \n", G.x2b(j),P(nden,kc,j,ic),P(nene,kc,j,ic) );
  }
  std::fclose(fp);

  // for z
  std::snprintf(fname, sizeof(fname), "snap/snz_%02d_%05d.dat",myid_w, index);
  fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int k=ks-2;k<=ke+2;k++){
    std::fprintf(fp, "%e %e %e \n", G.x3b(k),P(nden,k,jc,ic),P(nene,k,jc,ic) );
  }
  std::fclose(fp);
  
  index += 1;
}

void Output(bool& forcedamp){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  namespace io = mpi_dataio_mod;
  static int index = 0;
  
  static bool is_inited = false;

  
  if(!forcedamp && time_sim < time_out + dtout) return;
  
  if(myid_w==0) printf("output index=%i, time=%e \n",index,time_sim);
  
#pragma omp target update from (P.data[0:P.size])
  
  if (! is_inited){
    std::strcpy(io::id,"DT");
    std::strcpy(io::datadir,"./bindata/");
    (void)system("mkdir -p bindata");

    io::ntotal[dir1] = ntiles[dir1]*ngrid1;
    io::ntotal[dir2] = ntiles[dir2]*ngrid2;
    io::ntotal[dir3] = ntiles[dir3]*ngrid3;
    
    io::npart[dir1]  = ngrid1;
    io::npart[dir2]  = ngrid2;
    io::npart[dir3]  = ngrid3;
    
    io::nvars = 10;
    io::nvarg = 2;
    int ngrid1mod = ngrid1;
    int ngrid2mod = ngrid2;
    int ngrid3mod = ngrid3;
    // for grid, we include the edge
    if(coords[dir1] == ntiles[dir1] -1) ngrid1mod += 1;
    if(coords[dir2] == ntiles[dir2] -1) ngrid2mod += 1;
    if(coords[dir3] == ntiles[dir3] -1) ngrid3mod += 1;
    
    io::gridXout.allocate(io::nvarg,ngrid1mod);
    io::gridYout.allocate(io::nvarg,ngrid2mod);
    io::gridZout.allocate(io::nvarg,ngrid3mod);
    io::Fieldout.allocate(io::nvars,ngrid3,ngrid2,ngrid1);

    for(int i=0;i<ngrid1mod;i++){
      io::gridXout(0,i) = G.x1b(is+i);
      io::gridXout(1,i) = G.x1a(is+i);
    }
    
    for(int j=0;j<ngrid2mod;j++){
      io::gridYout(0,j) = G.x2b(js+j);
      io::gridYout(1,j) = G.x2a(js+j);
    }
    for(int k=0;k<ngrid3mod;k++){
      io::gridZout(0,k) = G.x3b(ks+k);
      io::gridZout(1,k) = G.x3a(ks+k);
    }
    
    is_inited = true;
  };
  
  // ---- output text (unf%05d.dat) ----
  char fname_unf[256];
  std::snprintf(fname_unf, sizeof(fname_unf), "bindata/unf%05d.dat", index);
  FILE* fp_unf = std::fopen(fname_unf, "w");
  if (!fp_unf){
    std::fprintf(stderr, "open failed: %s : %s\n", fname_unf, std::strerror(errno));
  }
  
  std::fprintf(fp_unf, "# %.17g %.17g\n", time_sim,dt);
  std::fprintf(fp_unf, "# %d \n", ngrid1);
  std::fprintf(fp_unf, "# %d \n", ngrid2);
  std::fprintf(fp_unf, "# %d \n", ngrid3);
  std::fclose(fp_unf);
  
  // ---- output data (bin%05d.dat) ----

  for (int k=0;k<ngrid3;k++)
    for (int j=0;j<ngrid2;j++)
      for (int i=0;i<ngrid1;i++){
	io::Fieldout(0,k,j,i) = P(nden,k+ks,j+js,i+is);
	io::Fieldout(1,k,j,i) = P(nve1,k+ks,j+js,i+is);
	io::Fieldout(2,k,j,i) = P(nve2,k+ks,j+js,i+is);
	io::Fieldout(3,k,j,i) = P(nve3,k+ks,j+js,i+is);
	io::Fieldout(4,k,j,i) = P(nbm1,k+ks,j+js,i+is);
	io::Fieldout(5,k,j,i) = P(nbm2,k+ks,j+js,i+is);
	io::Fieldout(6,k,j,i) = P(nbm3,k+ks,j+js,i+is);
	io::Fieldout(7,k,j,i) = P(nbps,k+ks,j+js,i+is);
	io::Fieldout(8,k,j,i) = P(npre,k+ks,j+js,i+is);
	io::Fieldout(9,k,j,i) = P(nxc ,k+ks,j+js,i+is);
  }
  io::MPIOutputBindary(index);
  WriteXMF(index, time_sim);
  index += 1;
}


int main() {
  
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace boundary_mod;
  using namespace mpi_config_mod;
  const bool NoOutput = false;
  static bool is_final = false;

  periodic[dir1] = 1;
  periodic[dir2] = 1;
  periodic[dir3] = 1;
  ntiles[dir1]   = 1;
  ntiles[dir2]   = 2;
  ntiles[dir3]   = 1;
  InitializeMPI();
  
  if(myid_w == 0) printf("setup grids and fields\n");
  
  AllocateHydroVariables(G,U,Fx,Fy,Fz,P);
 
  AllocateBoundaryVariables(Bs,Br);
  
  if (myid_w == 0) printf("grid size for x y z = %i %i %i\n",ngrid1*ntiles[dir1],ngrid2*ntiles[dir2],ngrid3*ntiles[dir3]);
  
  GenerateGrid(G);
  GenerateProblem(G,P,U);
  if (myid_w == 0) printf("entering main loop\n");
  int step = 0;
  auto time_beg = std::chrono::high_resolution_clock::now();

  for (step=0;step<stepmax;step++){
    ControlTimestep(G); 
    if (myid_w==0 && step%300 ==0 && !NoOutput) printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
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
    if (! NoOutput) Output(is_final);
    //if (! NoOutput) Output1D(is_final);

    if(time_sim > time_max) break;
    
  }

  //DeallocateHydroVariables(U,Fx,Fy,Fz,P);
  //DeallocateBoundaryVariables(Xs,Xe,Ys,Ye,Zs,Ze);
  
  auto time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = time_end - time_beg;
  if (myid_w == 0) printf("exiting main loop time=%e, step=%i\n",time_sim,step);
  if (myid_w == 0) printf("sim time [s]: %e\n", elapsed.count());
  if (myid_w == 0) printf("time/count/cell : %e\n", elapsed.count()/(ngrid1*ngrid2*ngrid3)/stepmax);

  is_final = true;
  //Output(is_final);
  //Output1D(is_final);
  
  printf("program has been finished\n");
  return 0;
}
