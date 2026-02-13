
#include<string>
#include<cstdio>
#include<cstring>
#include<cerrno>
#include<cstdint>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<vector>
#include<sys/stat.h>
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

  // --------------------------------------------------
  // small utilities (Fortran makedirs/to_str equivalents)
  // --------------------------------------------------
  static inline void makedirs(const std::string& dir){
    // portable enough for typical HPC linux environments
    (void) ::mkdir(dir.c_str(), 0755);
  }

  static inline std::string to_str_int(int v){
    std::ostringstream os;
    os << v;
    return os.str();
  }

  static inline std::string zpad_int(int v, int width){
    std::ostringstream os;
    os << std::setw(width) << std::setfill('0') << v;
    return os.str();
  }

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

  // ==================================================
  // WriteXDMF (Fortran: subroutine WriteXDMF)
  // ==================================================
  void WriteXDMF(double time, int nout){
    using namespace mpi_config_mod;

    // only rank0 writes
    if(myid_w != 0) return;

    makedirs("bindata");

    const int itot = ntotal[dir1];
    const int jtot = ntotal[dir2];
    const int ktot = ntotal[dir3];

    const std::string dirname = std::string(datadir);
    const std::string xmfname = dirname + "field" + zpad_int(nout,5) + ".xmf";
    const std::string fgridx  = "grid1D.bin";
    const std::string fgridy  = "grid2D.bin";
    const std::string fgridz  = "grid3D.bin";
    const std::string fdata   = "field" + zpad_int(nout,5) + ".bin";

    const std::int64_t bytes_per_real  = static_cast<std::int64_t>(sizeof(double));
    const std::int64_t ncell           = static_cast<std::int64_t>(itot) * jtot * ktot;
    const std::int64_t bytes_per_field = ncell * bytes_per_real;

    std::ofstream u(xmfname, std::ios::out | std::ios::trunc);
    if(!u){
      std::fprintf(stderr, "WriteXDMF: failed to open %s\n", xmfname.c_str());
      return;
    }

    u << "<?xml version=\"1.0\" ?>\n";
    u << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    u << "<Xdmf Version=\"2.0\">\n";
    u << "  <Domain>\n";
    u << "    <Grid Name=\"Grid\" GridType=\"Uniform\">\n";
    u << "      <Time Value=\"" << std::setprecision(16) << std::scientific << time << "\"/>\n";

    // NOTE: Our MPI-IO layout is C-order with global dims [nvars][Nx][Ny][Nz].
    // To keep the XDMF consistent with the binary ordering (last dimension fastest),
    // we use topology/attributes ordering as (Nx,Ny,Nz).
    u << "      <Topology TopologyType=\"3DRectMesh\" Dimensions=\""
      << (itot+1) << " " << (jtot+1) << " " << (ktot+1) << "\"/>\n";

    u << "      <Geometry GeometryType=\"VXVYVZ\">\n";
    // grid*.bin stores two 1D arrays: (cell edge, cell center). we reference the 2nd column.
    auto write_axis = [&](const std::string& fname, int n){
      const std::int64_t seek = static_cast<std::int64_t>(n) * bytes_per_real; // skip x?b
      u << "        <DataItem Dimensions=\"" << n << "\" NumberType=\"Float\" Precision=\"8\" "
           "Format=\"Binary\" Endian=\"Little\" Seek=\"" << seek << "\">"
        << fname << "</DataItem>\n";
    };
    write_axis(fgridx, itot+1);
    write_axis(fgridy, jtot+1);
    write_axis(fgridz, ktot+1);
    u << "      </Geometry>\n";

    auto write_attr = [&](const std::string& name, std::int64_t seek){
      u << "      <Attribute Name=\"" << name << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
      u << "        <DataItem Dimensions=\"" << itot << " " << jtot << " " << ktot
        << "\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\" Seek=\""
        << seek << "\">" << fdata << "</DataItem>\n";
      u << "      </Attribute>\n";
    };

    // variable names consistent with MPI_IO_Pack order.
    // (C++ side may later extend to gp and Xcomp; this keeps a Fortran-like convention.)
    static const char* base_names[] = {"d","v1","v2","v3","b1","b2","b3","bp","p"};
    std::int64_t off = 0;
    for(int n=0; n<nvars; ++n){
      std::string vname;
      if(n < 9) vname = base_names[n];
      else      vname = "X" + to_str_int(n-8); // X1, X2, ... (Fortran uses X1..Xncomp)
      write_attr(vname, off);
      off += bytes_per_field;
    }

    u << "    </Grid>\n";
    u << "  </Domain>\n";
    u << "</Xdmf>\n";
  }

  // ==================================================
  // ASC_Write (Fortran: subroutine ASC_WRITE)
  // ==================================================
  void ASC_Write(int nout){
    using namespace mpi_config_mod;
    using namespace resolution_mod;
    using namespace hydflux_mod;

    static bool is_inited_asc = false;
    if(!is_inited_asc){
      makedirs("ascdata");
      is_inited_asc = true;
    }

    const std::string fname = std::string("ascdata/") + "snap" + zpad_int(myid_w,3) + "-" + zpad_int(nout,5) + ".csv";
    std::ofstream f(fname, std::ios::out | std::ios::trunc);
    if(!f){
      std::fprintf(stderr, "ASC_Write: failed to open %s\n", fname.c_str());
      return;
    }

    const int k = ks;
    f << "# " << std::setprecision(6) << std::scientific << time_sim << " " << dt << "\n";
    f << "# " << (ie-is+1) << "\n";
    f << "# " << (je-js+1) << "\n";
    f << "# " << k << " " << std::setprecision(6) << std::scientific << G.x3b(k) << "\n";

    // Fortran header: "# x y d vx vy p phi" + Xcomp...
    // C++版は phi/gp をまだ持たないので、現状ある変数だけを出力。
    f << "# x y d vx vy p";
    // If extra scalars exist later, keep Fortran-like naming.
    for(int n=9; n<nvars; ++n) f << " X" << (n-8);
    f << "\n";

    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
        f << std::setprecision(6) << std::scientific
          << G.x1b(i) << ' ' << G.x2b(j) << ' '
          << P(nden,k,j,i) << ' ' << P(nve1,k,j,i) << ' ' << P(nve2,k,j,i) << ' ' << P(npre,k,j,i);

        // Append extra scalar(s) if they exist in Fieldout (packed) and we keep same order.
        // Note: in current code nvars==9 so this loop is a no-op.
        for(int n=9; n<nvars; ++n){
          // No direct primitive slot for Xcomp in this snapshot format yet.
          // Users can extend by adding the primitive index and printing it here.
          f << ' ' << 0.0;
        }
        f << "\n";
      }
      f << "\n";
    }
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

  // ==================================================
  // WriteXDMF / ASC_Write  (ported from output.f90)
  // ==================================================

  namespace {
    inline void mkdir_p(const char* dir) {
      // best-effort: ignore if already exists
      (void)mkdir(dir, 0777);
    }

    inline std::string zpad_int(int v, int width) {
      std::ostringstream ss;
      ss << std::setw(width) << std::setfill('0') << v;
      return ss.str();
    }

    inline void write_axis(std::ofstream& ofs,
                           const std::string& fname,
                           int n,
                           std::int64_t seek_bytes,
                           int bytes_per_real) {
      ofs << "        <DataItem Dimensions=\"" << n
          << "\" NumberType=\"Float\" Precision=\"" << bytes_per_real
          << "\" Format=\"Binary\" Endian=\"Little\" Seek=\"" << seek_bytes
          << "\"  >" << fname << "</DataItem>\n";
    }

    inline void write_attr(std::ofstream& ofs,
                           const std::string& name,
                           const std::string& fname,
                           int nx, int ny, int nz,
                           std::int64_t seek_bytes,
                           int bytes_per_real) {
      ofs << "      <Attribute Name=\"" << name
          << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
      // NOTE: our MPI-IO writes a global array of shape [nvars][nx][ny][nz] in C-order.
      // For XDMF readers that assume the last dimension is fastest, we keep Dimensions="nx ny nz".
      ofs << "        <DataItem Dimensions=\"" << nx << " " << ny << " " << nz
          << "\" NumberType=\"Float\" Precision=\"" << bytes_per_real
          << "\" Format=\"Binary\" Endian=\"Little\" Seek=\"" << seek_bytes
          << "\" >" << fname << "</DataItem>\n";
      ofs << "      </Attribute>\n";
    }

    inline std::vector<std::string> default_varnames(int nvars_) {
      // Matches the order used in MPI_IO_Pack
      std::vector<std::string> names;
      names.reserve(static_cast<size_t>(nvars_));
      names.push_back("d");
      names.push_back("v1");
      names.push_back("v2");
      names.push_back("v3");
      names.push_back("b1");
      names.push_back("b2");
      names.push_back("b3");
      names.push_back("bp");
      names.push_back("p");
      // Any remaining variables are treated as passive scalars (X1, X2, ...)
      for (int n = static_cast<int>(names.size()); n < nvars_; ++n) {
        const int idx = n - 8; // so first extra becomes X1
        names.push_back("X" + std::to_string(idx));
      }
      if ((int)names.size() > nvars_) names.resize(nvars_);
      return names;
    }
  } // anonymous namespace

  void WriteXDMF(double time, int nout) {
    using namespace resolution_mod;

    mkdir_p(datadir);

    const int itot = ntotal[0];
    const int jtot = ntotal[1];
    const int ktot = ntotal[2];

    const int bytes_per_real = (int)sizeof(double);
    const std::int64_t ncell = (std::int64_t)itot * (std::int64_t)jtot * (std::int64_t)ktot;
    const std::int64_t bytes_per_field = ncell * (std::int64_t)bytes_per_real;

    const std::string xmfname = std::string(datadir) + "field" + zpad_int(nout,5) + ".xmf";
    const std::string fgridx  = "grid1D.bin";
    const std::string fgridy  = "grid2D.bin";
    const std::string fgridz  = "grid3D.bin";
    const std::string fdata   = "field" + zpad_int(nout,5) + ".bin";

    std::ofstream ofs(xmfname);
    if (!ofs) {
      std::fprintf(stderr, "WriteXDMF: open failed: %s\n", xmfname.c_str());
      return;
    }

    ofs << "<?xml version=\"1.0\" ?>\n";
    ofs << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    ofs << "<Xdmf Version=\"2.0\">\n";
    ofs << "  <Domain>\n";
    ofs << "    <Grid Name=\"Grid\" GridType=\"Uniform\">\n";
    ofs << "      <Time Value=\"" << std::setprecision(16) << time << "\"/>\n";

    // NOTE:
    // Our global array ordering for fields is [nvars][nx][ny][nz] (C-order).
    // We therefore set Topology/Attribute Dimensions in (nx,ny,nz) order.
    ofs << "      <Topology TopologyType=\"3DRectMesh\" Dimensions=\""
        << (itot + 1) << " " << (jtot + 1) << " " << (ktot + 1)
        << "\"/>\n";

    ofs << "      <Geometry GeometryType=\"VXVYVZ\">\n";
    // grid files store [x_b][x_a] (2 arrays). We want x_a (cell edges) => Seek = (itot+1)*8
    write_axis(ofs, fgridx, itot + 1, (std::int64_t)(itot + 1) * bytes_per_real, bytes_per_real);
    write_axis(ofs, fgridy, jtot + 1, (std::int64_t)(jtot + 1) * bytes_per_real, bytes_per_real);
    write_axis(ofs, fgridz, ktot + 1, (std::int64_t)(ktot + 1) * bytes_per_real, bytes_per_real);
    ofs << "      </Geometry>\n";

    std::int64_t off_base = 0;
    const auto names = default_varnames(nvars);
    for (int n = 0; n < nvars; ++n) {
      write_attr(ofs, names[(size_t)n], fdata, itot, jtot, ktot, off_base, bytes_per_real);
      off_base += bytes_per_field;
    }

    ofs << "    </Grid>\n";
    ofs << "  </Domain>\n";
    ofs << "</Xdmf>\n";
  }

  void ASC_Write(int nout) {
    using namespace mpi_config_mod;
    using namespace resolution_mod;
    using namespace hydflux_mod;

    static bool is_inited_local = false;
    if (!is_inited_local) {
      mkdir_p("ascdata");
      is_inited_local = true;
    }

    std::ostringstream fn;
    fn << "ascdata/snap" << std::setw(3) << std::setfill('0') << myid_w
       << "-" << std::setw(5) << std::setfill('0') << nout << ".csv";

    std::ofstream ofs(fn.str());
    if (!ofs) {
      std::fprintf(stderr, "ASC_Write: open failed: %s\n", fn.str().c_str());
      return;
    }

    const int k = ks;
    ofs << "# " << std::scientific << std::setprecision(6) << time_sim << " " << dt << "\n";
    ofs << "# " << (ie - is + 1) << "\n";
    ofs << "# " << (je - js + 1) << "\n";
    ofs << "# " << k << " " << std::scientific << std::setprecision(6) << G.x3b(k) << "\n";
    ofs << "# x y d vx vy p\n";

    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        ofs << std::scientific << std::setprecision(6)
            << G.x1b(i) << " "
            << G.x2b(j) << " "
            << P(nden, k, j, i) << " "
            << P(nve1, k, j, i) << " "
            << P(nve2, k, j, i) << " "
            << P(npre, k, j, i)
            << "\n";
      }
      ofs << "\n";
    }
  }

}// namespace mpi_dataio_mod

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
    if(myid_w==0) mpiio::WriteXDMF(time_sim,index);
  }else{
    mpiio::ASC_Write(index);
  };
  if(myid_w==0) printf("output index=%i, time=%e \n",index,time_sim);
  
  index += 1;
  time_out = time_sim;
}
