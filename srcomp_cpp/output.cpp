
#include<string>
#include<cstdio>
#include<cstring>
#include<cerrno>
#include "mpi_config.hpp"
#include "output.hpp"
#include "hydro.hpp"

#include "resolution.hpp"

// ==================================================
// DATA IO (MPI-IO)
// ==================================================
namespace mpi_dataio_mod {

  int ntotal[3] = {0,0,0};
  int npart[3]  = {0,0,0};
  int nvars = 0, nvarg = 0;

  // Like Fortran (output.f90):
  //   logical,parameter :: binaryout = .true.
  // Users can switch this from main.
  bool binaryout = true;

  GridArray<double> gridXout;
  GridArray<double> gridYout;
  GridArray<double> gridZout;

  FieldArray<double> Fieldout;
  
  static MPI_Datatype SAG1D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAG2D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAG3D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAD3D = MPI_DATATYPE_NULL;
  
  static bool is_inited = false;

  char id[10];
  char datadir[10];

  // --------------------------------------------------
  // ASCII output (Fortran: ASC_WRITE)
  //   - quick check / header-only output
  //   - writes `bindata/unf%05d.dat` (one file per rank)
  //     *# time dt
  //     *# ngrid1_local
  //     *# ngrid2_local
  //     *# ngrid3_local
  // --------------------------------------------------
  void ASC_WRITE(int timeid, double time, double dt,
                 int ngrid1_local, int ngrid2_local, int ngrid3_local) {
    // Ensure directory exists (same spirit as Fortran's makedirs)
    (void)system("mkdir -p bindata");

    char fname[256];
    std::snprintf(fname, sizeof(fname), "bindata/unf%05d.dat", timeid);
    FILE* fp = std::fopen(fname, "w");
    if (!fp) {
      std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
      return;
    }
    std::fprintf(fp, "# %.17g %.17g\n", time, dt);
    std::fprintf(fp, "# %d \n", ngrid1_local);
    std::fprintf(fp, "# %d \n", ngrid2_local);
    std::fprintf(fp, "# %d \n", ngrid3_local);
    std::fclose(fp);
  }

// ---- MPI-IO 書き出し本体（Fortran: MPIOutputBindary）----
  void MPIOutputBindary(int timeid) {
  using namespace mpi_config_mod;

  int Asize[4], Ssize[4], Start[4];
  MPI_Offset idisp = 0;

  // ============ Grid 1D ============
  if (!is_inited) {
    Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
    Asize[1] = ntotal[dir1] + 1;
	    Ssize[1] = npart[dir1];
    if(coords[dir1]==ntiles[dir1]-1) Ssize[1] += 1;
    Start[1] = npart[dir1] * coords[dir1];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_C,MPI_DOUBLE, &SAG1D);
    MPI_Type_commit(&SAG1D);

    std::string fpath = std::string(datadir) + "g1d" + id;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG1D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count = static_cast<MPI_Offset>(Ssize[1]) * nvarg;
    MPI_File_write_all(fh, gridXout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // ============ Grid 2D ============
  if (!is_inited) {
    Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
    Asize[1] = ntotal[dir2] + 1;
	    Ssize[1] = npart[dir2];
    if(coords[dir2]==ntiles[dir2]-1) Ssize[1] += 1;
	    Start[1] = npart[dir2] * coords[dir2];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_C,MPI_DOUBLE, &SAG2D);
    MPI_Type_commit(&SAG2D);

    std::string fpath = std::string(datadir) + "g2d" + id;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG2D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count = static_cast<MPI_Offset>(Ssize[1]) * nvarg;
    MPI_File_write_all(fh, gridYout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // ============ Grid 3D ============
  if (!is_inited) {
    Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
    Asize[1] = ntotal[dir3] + 1;
	    Ssize[1] = npart[dir3];
    if(coords[dir3]==ntiles[dir3]-1) Ssize[1] += 1;
    Start[1] = npart[dir3] * coords[dir3];

	    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_C, MPI_DOUBLE, &SAG3D);
    MPI_Type_commit(&SAG3D);

    std::string fpath = std::string(datadir) + "g3d" + id;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG3D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count = static_cast<MPI_Offset>(Ssize[1]) * nvarg;
    MPI_File_write_all(fh, gridZout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // ============ Data 3D タイプ準備 ============
  if (!is_inited) {
    Asize[0] = nvars; Ssize[0] = nvars; Start[0] = 0;
    Asize[1] = ntotal[0]; Ssize[1] = npart[0]; Start[1] = npart[0] * coords[0];
    Asize[2] = ntotal[1]; Ssize[2] = npart[1]; Start[2] = npart[1] * coords[1];
    Asize[3] = ntotal[2]; Ssize[3] = npart[2]; Start[3] = npart[2] * coords[2];

	    MPI_Type_create_subarray(4, Asize, Ssize, Start, MPI_ORDER_C,
                             MPI_DOUBLE, &SAD3D);
    MPI_Type_commit(&SAD3D);
  }

  // ============ Data 書き出し ============
  {
    char filename[64];
    std::snprintf(filename, sizeof(filename), "d3d%s.%05d", id, timeid);
    std::string fpath = std::string(datadir) + filename;

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAD3D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count =
      static_cast<MPI_Offset>(nvars) *
      static_cast<MPI_Offset>(npart[0]) *
      static_cast<MPI_Offset>(npart[1]) *
      static_cast<MPI_Offset>(npart[2]);

    MPI_File_write_all(fh, Fieldout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  is_inited = true;
}

} // namespace mpiiomod

void Output1D(bool is_final){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  static int index = 0;  
  static bool is_inited = false;
  const int dir1=0, dir2=1, dir3=2;
  // Follow Fortran's output.f90:
  //   if(time < tout+dtout .and. .not. is_final) return
  if(!is_final && time_sim < time_out + dtout) return;

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
  // Record last output time (Fortran: tout=time)
  time_out = time_sim;
}

void Output(bool is_final){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  namespace io = mpi_dataio_mod;
  static int index = 0;
  
  static bool is_inited = false;

  
  // Follow Fortran's output.f90:
  //   if(time < tout+dtout .and. .not. is_final) return
  if(!is_final && time_sim < time_out + dtout) return;
  
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
    
    io::nvars = 9;
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

    for(int i=0;i<=ngrid1mod;i++){
      io::gridXout(0,i) = G.x1b(is+i);
      io::gridXout(1,i) = G.x1a(is+i);
    }
    
    for(int j=0;j<=ngrid2mod;j++){
      io::gridYout(0,j) = G.x2b(js+j);
      io::gridYout(1,j) = G.x2a(js+j);
    }
    for(int k=0;k<=ngrid3mod;k++){
      io::gridZout(0,k) = G.x3b(ks+k);
      io::gridZout(1,k) = G.x3a(ks+k);
    }
    
    is_inited = true;
  };
  
  // Follow Fortran's Output() logic:
  //   if(binaryout) then
  //     MPI binary output
  //   else
  //     ASCII output
  //   endif
  if(!io::binaryout){
    io::ASC_WRITE(index, time_sim, dt, ngrid1, ngrid2, ngrid3);
    index += 1;
    time_out = time_sim;
    return;
  }

  // ---- output data (MPI binary) ----

  for (int k=0;k<=ngrid3;k++)
    for (int j=0;j<=ngrid2;j++)
      for (int i=0;i<=ngrid1;i++){
	io::Fieldout(0,k,j,i) = P(nden,k+ks,j+js,i+is);
	io::Fieldout(1,k,j,i) = P(nve1,k+ks,j+js,i+is);
	io::Fieldout(2,k,j,i) = P(nve2,k+ks,j+js,i+is);
	io::Fieldout(3,k,j,i) = P(nve3,k+ks,j+js,i+is);
	io::Fieldout(4,k,j,i) = P(nbm1,k+ks,j+js,i+is);
	io::Fieldout(5,k,j,i) = P(nbm2,k+ks,j+js,i+is);
	io::Fieldout(6,k,j,i) = P(nbm3,k+ks,j+js,i+is);
	io::Fieldout(7,k,j,i) = P(nbps,k+ks,j+js,i+is);
	io::Fieldout(8,k,j,i) = P(npre,k+ks,j+js,i+is); 
  }
  io::MPIOutputBindary(index);
  index += 1;
  time_out = time_sim;
}
