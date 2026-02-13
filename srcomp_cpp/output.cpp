
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

  GridArray<double> gridXout;
  GridArray<double> gridYout;
  GridArray<double> gridZout;

  FieldArray<double> Fieldout;
  
  static MPI_Datatype SAG1D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAG2D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAG3D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAD3D = MPI_DATATYPE_NULL;
  
  static bool is_inited = false;

  bool binaryout = true;
  
  char datadir[10] = "bindata/"; 

  void MPI_IO_Write(int timeid) {
    using namespace mpi_config_mod;
    using namespace resolution_mod;
    int Asize[4], Ssize[4], Start[4];
    MPI_Offset idisp = 0;

    char fname[256];
    std::snprintf(fname, sizeof(fname), "bindata/unf%05d.dat", timeid);
    FILE* fp = std::fopen(fname, "w");
    if (!fp) {
      std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
      return;
    }
    std::fprintf(fp, "# %.17g %.17g\n", time, dt);
    std::fprintf(fp, "# %d \n", ntotal[dir1]);
    std::fprintf(fp, "# %d \n", ntotal[dir2]);
    std::fprintf(fp, "# %d \n", ntotal[dir3]);
    std::fclose(fp);
    
    
    // ============ Grid 1D ============
    if (!is_inited) {
      Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
      Asize[1] = ntotal[dir1] + 1;
      Ssize[1] = npart[dir1];
      if(coords[dir1]==ntiles[dir1]-1) Ssize[1] += 1;
      Start[1] = npart[dir1] * coords[dir1];

      MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_C,MPI_DOUBLE, &SAG1D);
      MPI_Type_commit(&SAG1D);

      std::string fpath = std::string(datadir) + "grid1D.bin";
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

      std::string fpath = std::string(datadir) + "grid2d.bin";
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

      std::string fpath = std::string(datadir) + "grid3D.bin";
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
    std::snprintf(filename, sizeof(filename), "field%05d.bin", timeid);
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

  void MPI_IO_Pack(int index){
    using namespace mpi_config_mod;
    using namespace resolution_mod;
    using namespace hydflux_mod;
    
    static bool is_inited = false;
  
    if (! is_inited){
      (void)system("mkdir -p bindata");


      ntotal[dir1] = ntiles[dir1]*ngrid1;
      ntotal[dir2] = ntiles[dir2]*ngrid2;
      ntotal[dir3] = ntiles[dir3]*ngrid3;
    
      npart[dir1]  = ngrid1;
      npart[dir2]  = ngrid2;
      npart[dir3]  = ngrid3;
    
      nvars = 9;
      nvarg = 2;
      int ngrid1mod = ngrid1;
      int ngrid2mod = ngrid2;
      int ngrid3mod = ngrid3;
      // for grid, we include the edge
      if(coords[dir1] == ntiles[dir1] -1) ngrid1mod += 1;
      if(coords[dir2] == ntiles[dir2] -1) ngrid2mod += 1;
      if(coords[dir3] == ntiles[dir3] -1) ngrid3mod += 1;
    
      gridXout.allocate(nvarg,ngrid1mod);
      gridYout.allocate(nvarg,ngrid2mod);
      gridZout.allocate(nvarg,ngrid3mod);
      Fieldout.allocate(nvars,ngrid3,ngrid2,ngrid1);

      for(int i=0;i<=ngrid1mod;i++){
	gridXout(0,i) = G.x1b(is+i);
	gridXout(1,i) = G.x1a(is+i);
      }
    
      for(int j=0;j<=ngrid2mod;j++){
	gridYout(0,j) = G.x2b(js+j);
	gridYout(1,j) = G.x2a(js+j);
      }
      for(int k=0;k<=ngrid3mod;k++){
	gridZout(0,k) = G.x3b(ks+k);
	gridZout(1,k) = G.x3a(ks+k);
      }
    
      is_inited = true;
    };

  // ---- output data (MPI binary) ----

    for (int k=0;k<=ngrid3;k++)
    for (int j=0;j<=ngrid2;j++)
    for (int i=0;i<=ngrid1;i++){
      Fieldout(0,k,j,i) = P(nden,k+ks,j+js,i+is);
      Fieldout(1,k,j,i) = P(nve1,k+ks,j+js,i+is);
      Fieldout(2,k,j,i) = P(nve2,k+ks,j+js,i+is);
      Fieldout(3,k,j,i) = P(nve3,k+ks,j+js,i+is);
      Fieldout(4,k,j,i) = P(nbm1,k+ks,j+js,i+is);
      Fieldout(5,k,j,i) = P(nbm2,k+ks,j+js,i+is);
      Fieldout(6,k,j,i) = P(nbm3,k+ks,j+js,i+is);
      Fieldout(7,k,j,i) = P(nbps,k+ks,j+js,i+is);
      Fieldout(8,k,j,i) = P(npre,k+ks,j+js,i+is); 
  };

  }// MPI_IO_PACK
}// namespace mpiiomod

void Output(bool is_forced){
  using namespace mpi_config_mod;
  using namespace resolution_mod;
  using namespace hydflux_mod;
  namespace mpiio = mpi_dataio_mod;
  static int index = 0;
    
  if(!is_forced && time_sim < time_out + dtout) return;
  
  
#pragma omp target update from (P.data[0:P.size])

  //   if(binaryout) then
  //     MPI binary output
  //   else
  //     ASCII output
  //   endif
  if(mpiio::binaryout){
    mpiio::MPI_IO_Pack(index);
    mpiio::MPI_IO_Write(index);
    //    if(myid_w==0) mpiio:WriteXDMF(time,index);
  }else{
    // ASC_WRITE(index);
  };
  if(myid_w==0) printf("output index=%i, time=%e \n",index,time_sim);
  
  index += 1;
  time_out = time_sim;
}
