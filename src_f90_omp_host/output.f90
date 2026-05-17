
!==================================================
! DATA IO
!==================================================

module mpiiomod
  use basicmod, only: ncomp
  implicit none
  private
  integer::SAG1D,SAG2D,SAG3D,SAD3D
  integer :: commG1D,commG2D,commG3D
  
  integer,dimension(3):: ntotal
  integer,dimension(3):: npart
  integer,parameter:: nvarg = 2 !! number of the variables for grid
  integer,parameter:: nvars = 10+ncomp!! number of the variables for variable

  character(len=10),parameter :: binaryextension=".bin"
  character(len=10),parameter :: filegrid1D="grid1D"//binaryextension
  character(len=10),parameter :: filegrid2D="grid2D"//binaryextension
  character(len=10),parameter :: filegrid3D="grid3D"//binaryextension
  character(len=10),parameter :: filefield ="field"
  
  character(len=10),parameter :: datadir="bindata/"
    
  real(8),dimension(:,:),allocatable:: gridX, gridY, gridZ
  real(8),dimension(:,:,:,:),allocatable:: data3D
  public:: ntotal,npart
  public:: nvars,nvarg
  
  public:: gridX,gridY,gridZ,data3D
  public:: MPI_IO_Pack, MPI_IO_Write, WriteXDMF
contains  
  subroutine MPI_IO_Write(timeid)
    use mpimod
    implicit NONE
    integer,intent(in) :: timeid
    integer :: i, j, k, l, m, n

    integer::iss,jss,kss
    integer::iee,jee,kee
    integer::itot,jtot,ktot
      
    integer strtoi 
    character(len=40)::usrfile
    character(len=30) :: fpathbin,fpathunf
    integer,save:: unitd3d,unitg1d, unitg2d, unitg3d, unitg0d
    data unitd3d / 512 /
    data unitg1d / 513 /
    data unitg2d / 514 /
    data unitg3d / 515 /
    data unitg0d / 516 /

    logical :: fileflag

    logical,save :: is_inited
    data is_inited / .false. /
    
    integer,dimension(4)::Asize,Ssize,Start
    integer(kind=MPI_OFFSET_KIND) idisp
    data idisp / 0 /

    integer::color,key
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    init1D: if(.not. is_inited )then
       Asize(2) = nvarg
       Ssize(2) = nvarg
       Start(2) = 0
       
       Asize(1) = ntotal(1)+1 ! total izones + edge
       Ssize(1) = npart(1) ! izones in 1 process
       if(coords(1) .eq. ntiles(1)-1)Ssize(1)=Ssize(1)+1  ! + edge
       Start(1) = npart(1) * coords(1)
              
       call MPI_TYPE_CREATE_SUBARRAY( &
     & 2, & ! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION,& 
     & SAG1D, &! Data type of Subarray for Grid 1D
     & ierr)
       
       call MPI_TYPE_COMMIT(SAG1D,ierr)

      color =  coords(2)*ntiles(3)+coords(3)
      key   =  coords(1)
      call MPI_COMM_SPLIT(comm3d,color,key,commG1D,ierr)
      
      color1D: if(color == 0) then
      write(usrfile,"(A)") trim(filegrid1D)
      fpathbin = trim(datadir)//usrfile
      call MPI_FILE_OPEN(commG1D, &
     &                         fpathbin, &  ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg1d,ierr)

      call MPI_FILE_SET_VIEW( &
     &   unitg1d, &! file path
     &     idisp, &! 
     & MPI_DOUBLE_PRECISION, & 
     &     SAG1D, &! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &   unitg1d, &  ! file path
     &     gridX, &  ! the data
     & Ssize(2)*Ssize(1), &! total data number
     & MPI_DOUBLE_PRECISION, & 
     & mpi_status_ignore, &
     & ierr)
      call MPI_FILE_CLOSE(unitg1d,ierr)
      endif color1D
   
   endif init1D
!   print *, "p2"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   init2D: if(.not. is_inited )then
      Asize(2) = nvarg
      Ssize(2) = nvarg
      Start(2) = 0
      
      Asize(1) = ntotal(2)+1 ! total jzones + edge
      Ssize(1) = npart(2)    ! jzones in 1 process
      if(coords(2) .eq. ntiles(2)-1)Ssize(1)=Ssize(1)+1  ! + edge
      Start(1) = npart(2) * coords(2)
      call MPI_TYPE_CREATE_SUBARRAY(&
     & 2, & ! dimension of array
     & Asize,Ssize,Start,&
     & MPI_ORDER_FORTRAN,&
     & MPI_DOUBLE_PRECISION,&
     & SAG2D, &! Data type of Subarray for Grid 2D
     & ierr)
      call MPI_TYPE_COMMIT(SAG2D,ierr)
         
      color =  coords(3)*ntiles(1)+coords(1)
      key   =  coords(2)
      call MPI_COMM_SPLIT(comm3d,color,key,commG2D,ierr)
      
      color2D: if(color == 0) then
      write(usrfile,"(A)") trim(filegrid2D)
      fpathbin = trim(datadir)//usrfile
 
      call MPI_FILE_OPEN(commG2D, &
     &                         fpathbin, & ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg2d,ierr)
      call MPI_FILE_SET_VIEW(&
     &   unitg2d, & ! file ID
     &     idisp, & ! 
     & MPI_DOUBLE_PRECISION,&
     &     SAG2D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &   unitg2d,  &! file ID
     &     gridY,  &! the data
     & Ssize(2)*Ssize(1),& ! total data number
     & MPI_DOUBLE_PRECISION,&
     & mpi_status_ignore,&
     & ierr)
      call MPI_FILE_CLOSE(unitg2d,ierr)
     endif color2D
   endif init2D
      
   init3D: if(.not. is_inited )then
      Asize(1) = ntotal(3)+1  ! total kzones+edge
      Ssize(1) = npart(3) ! kzones in 1 process

      Asize(2) = nvarg
      Ssize(2) = nvarg
      Start(2) = 0
      
      if(coords(3) .eq. ntiles(3)-1)Ssize(1)=Ssize(1)+1  ! + edge
      Start(1) = npart(3) * coords(3)
      call MPI_TYPE_CREATE_SUBARRAY(&
     & 2, & ! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION,&
     & SAG3D, &! Data type of Subarray for Grid 3D
     & ierr)
         call MPI_TYPE_COMMIT(SAG3D,ierr)

      color =  coords(1)*ntiles(2)+coords(2)
      key   =  coords(3)
      call MPI_COMM_SPLIT(comm3d,color,key,commG3D,ierr)
      color3D: if (color == 0) then
      write(usrfile,"(A)") trim(filegrid3D)
      fpathbin = trim(datadir)//usrfile
      
      call MPI_FILE_OPEN(commG3D,&
     &                         fpathbin,&  ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE,&
     &            MPI_INFO_NULL,unitg3d,ierr)
      call MPI_FILE_SET_VIEW(&
     &  unitg3d, &  ! file path
     &    idisp, &  ! 
     & MPI_DOUBLE_PRECISION, &
     &     SAG3D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL(&
     &  unitg3d, & ! file path
     &    gridZ, & ! the data
     & Ssize(2)*Ssize(1),& ! total data number
     & MPI_DOUBLE_PRECISION, &
     & mpi_status_ignore, &
     & ierr)
      call MPI_FILE_CLOSE(unitg3d,ierr)
      endif color3D
   endif init3D
          
    initdata: if(.not. is_inited )then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Asize(1) = ntotal(1) ! total zones for 1D 
       Asize(2) = ntotal(2) ! total zones for 2D
       Asize(3) = ntotal(3) ! total zones for 3D
       Ssize(1) =  npart(1) ! partial zones in 1 process 
       Ssize(2) =  npart(2) ! partial zones in 1 process 
       Ssize(3) =  npart(3) ! partial zones in 1 process 
       Start(1) =  npart(1) * coords(1)
       Start(2) =  npart(2) * coords(2)
       Start(3) =  npart(3) * coords(3)
       
       Asize(4) = nvars
       Ssize(4) = nvars
       Start(4) = 0 

       call MPI_TYPE_CREATE_SUBARRAY(&
     & 4, &! dimension of array
     & Asize,Ssize,Start,&
     & MPI_ORDER_FORTRAN,&
     & MPI_DOUBLE_PRECISION,&
     & SAD3D,& ! Data type of Subarray for data 3D
     & ierr)
       call MPI_TYPE_COMMIT(SAD3D,ierr)

    endif initdata
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(usrfile,"(A,i5.5,A)") trim(filefield),timeid,trim(binaryextension)
    fpathbin = trim(datadir)//usrfile
    
      call MPI_FILE_OPEN(MPI_COMM_WORLD, &
     &                         fpathbin, & ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitd3d,ierr)
      call MPI_FILE_SET_VIEW(&
     &  unitd3d,  &! file path
     &     idisp, & ! 
     & MPI_DOUBLE_PRECISION,& 
     &     SAD3D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL(&
     &   unitd3d,  &! file path
     &    data3D,  &! the data
     & npart(1)*npart(2)*npart(3)*nvars,& ! total data number
     & MPI_DOUBLE_PRECISION,&  
     & mpi_status_ignore,&
     &      ierr)
      call MPI_FILE_CLOSE(unitd3d,ierr)
      
      is_inited = .true.

      return
  end subroutine MPI_IO_Write
  
  subroutine MPI_IO_Pack(nout)
    !! Output the grid and variables.
    !! The grid data contain the information of cell center and cell edge.
    !! The variable data contains that of  cell center.
    use basicmod
    use mpimod
    implicit none
    integer::i,j,k,n
    integer,intent(in):: nout
    integer::iee,jee,kee
    character(20),parameter::dirname=datadir
    character(40)::filename
    integer,parameter::unitout=17
    integer,parameter:: gs=1
     
    logical, save:: is_inited
    data is_inited /.false./

    iee = ie
    jee = je
    kee = ke
  
    !> Include the information of the cell edge!
    if(coords(1) .eq. ntiles(1)-1) iee = ie+1
    if(coords(2) .eq. ntiles(2)-1) jee = je+1
    if(coords(3) .eq. ntiles(3)-1) kee = ke+1
    
    if (.not. is_inited) then
       npart(1) = ngrid1
       npart(2) = ngrid2
       npart(3) = ngrid3
       
       ntotal(1) = ngrid1*ntiles(1)
       ntotal(2) = ngrid2*ntiles(2)
       ntotal(3) = ngrid3*ntiles(3)
     
       allocate(gridX(1:iee-is+1,nvarg))
       allocate(gridY(1:jee-js+1,nvarg))
       allocate(gridZ(1:kee-ks+1,nvarg))

       allocate(data3D(ngrid1,ngrid2,ngrid3,nvars))
         
       call makedirs(dirname)
       is_inited =.true.
      endif
      
      !> Information of meta data.
      if(myid_w == 0)then
         write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
         filename = trim(dirname)//trim(filename)
         
         open(unitout,file=filename,status='replace',form='formatted')
         write(unitout,'(a2,2(1x,E12.3))') "# ",time,dt
         write(unitout,'(a2,1x,i5)') "# ",ngrid1*ntiles(1)
         write(unitout,'(a2,1x,i5)') "# ",ngrid2*ntiles(2)
         write(unitout,'(a2,1x,i5)') "# ",ngrid3*ntiles(3)
         close(unitout)
      endif
      
      gridX(1:iee-is+1,1) = x1b(is:iee) !! the final grid point is not necessary but outputed.  
      gridX(1:iee-is+1,2) = x1a(is:iee) !! the final grid is necessary. 
      
      gridY(1:jee-js+1,1) = x2b(js:jee) !! the final grid point is not necessary but outputed.
      gridY(1:jee-js+1,2) = x2a(js:jee) !! the final grid is necessary.
      
      gridZ(1:kee-ks+1,1) = x3b(ks:kee) !! the final grid point is not necessary but outputed.
      gridZ(1:kee-ks+1,2) = x3a(ks:kee) !! the final grid is necessary.

      !> The cell center value 
      data3D(1:npart(1),1:npart(2),1:npart(3), 1) =  d(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 2) = v1(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 3) = v2(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 4) = v3(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 5) = b1(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 6) = b2(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 7) = b3(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 8) = bp(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3), 9) =  p(is:ie,js:je,ks:ke)
      data3D(1:npart(1),1:npart(2),1:npart(3),10) = gp(is:ie,js:je,ks:ke)
      do n=1,ncomp
         data3D(1:npart(1),1:npart(2),1:npart(3),nvars-ncomp+n) = Xcomp(n,is:ie,js:je,ks:ke)
      enddo
      return
      !         write(6,*) "bpf2",nflux2(mbps,i,j,k)
  end subroutine MPI_IO_Pack

  subroutine WriteXDMF(time,nout)
    use, intrinsic :: iso_fortran_env, only: int64,real64
    implicit none
    real(8),intent(in)::time
    integer, intent(in) :: nout
    integer:: itot, jtot, ktot
    integer:: n

    character(*), parameter :: dirname = datadir
    character(256) :: xmfname, fgridx, fgridy, fgridz, fdata
    integer :: u, ncell
    integer(int64) :: bytes_per_real, bytes_per_field
    integer(int64) :: off_base


    itot = ntotal(1)
    jtot = ntotal(2)
    ktot = ntotal(3)
    
    ! ---- file names (match your ReadData) ----
    write(xmfname,'(a,i5.5,a)') trim(filefield), nout, ".xmf"
    xmfname = trim(dirname)//trim(xmfname)

    fgridx = trim(filegrid1D)
    fgridy = trim(filegrid2D)
    fgridz = trim(filegrid3D)
    write(fdata,'(a,i5.5,a)') trim(filefield), nout,trim(binaryextension)

    ! ---- sizes & offsets ----
    ! stream/unformatted wrote raw reals; assume real64 (8 bytes) because iso_fortran_env real64
    bytes_per_real = int(storage_size(0.0_real64)/8, int64)
    
    ncell = itot*jtot*ktot
    bytes_per_field = int(ncell, int64) * bytes_per_real

    ! ---- write XDMF (XML) ----
    open(newunit=u, file=xmfname, status="replace", action="write", form="formatted")
    
    write(u,'(a)') '<?xml version="1.0" ?>'
    write(u,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(u,'(a)') '<Xdmf Version="2.0">'
    write(u,'(a)') '  <Domain>'
    write(u,'(a)') '    <Grid Name="Grid" GridType="Uniform">'
    write(u,'(a,es24.16,a)') '      <Time Value="', time, '"/>'

    ! Topology: Rectilinear mesh uses node dimensions = (in+1, jn+1, kn+1)
    ! Use Order="Fortran" consistently so Dimensions list matches (i,j,k).
    write(u,'(a,i0,1x,i0,1x,i0,a)') '      <Topology TopologyType="3DRectMesh" Dimensions="', &
         ktot+1, jtot+1, itot+1, '"/>'

    write(u,'(a)') '      <Geometry GeometryType="VXVYVZ">'
    call write_axis(u, fgridx, itot+1, int(itot+1,int64)*bytes_per_real, bytes_per_real)  ! x1a(:)
    call write_axis(u, fgridy, jtot+1, int(jtot+1,int64)*bytes_per_real, bytes_per_real)  ! x2a(:)
    call write_axis(u, fgridz, ktot+1, int(ktot+1,int64)*bytes_per_real, bytes_per_real)  ! x3a(:)

    write(u,'(a)') '      </Geometry>'

    off_base = 0_int64
    ! ---- Cell-centered attributes (itot,jtot,ktot) ----
    call write_attr(u, "d" , fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "v1", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "v2", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "v3", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "b1", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "b2", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "b3", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "bp", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "p" , fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    call write_attr(u, "gp", fdata, itot, jtot, ktot, off_base, bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    do n=1,ncomp
       call write_attr(u, "X"//trim(adjustl(to_str(n))) , fdata, itot, jtot, ktot, off_base , bytes_per_real); off_base=off_base + 1_int64 * bytes_per_field
    enddo
    write(u,'(a)') '    </Grid>'
    write(u,'(a)') '  </Domain>'
    write(u,'(a)') '</Xdmf>'
    
    close(u)

  contains

    subroutine write_axis(u, fname, n, seek_bytes, bpr)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer, intent(in) :: u, n
      character(*), intent(in) :: fname
      integer(int64), intent(in) :: seek_bytes, bpr

      ! XDMF Binary DataItem:
      ! - Format="Binary"
      ! - Seek="bytes" to skip within file
      ! - Precision="8" etc.
      write(u,'(a)') '        <DataItem Dimensions="'//trim(itoa(n))// &
           '" NumberType="Float" Precision="'//trim(itoa(int(bpr)))// &
           '" Format="Binary" Endian="Little" Seek="'//trim(i64toa(seek_bytes))// &
           '"  >'//trim(fname)//'</DataItem>'
    end subroutine write_axis

    subroutine write_attr(u, name, fname, nx, ny, nz, seek_bytes, bpr)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer, intent(in) :: u, nx, ny, nz
      character(*), intent(in) :: name, fname
      integer(int64), intent(in) :: seek_bytes, bpr

      write(u,'(a)') '      <Attribute Name="'//trim(name)//'" AttributeType="Scalar" Center="Cell">'
      write(u,'(a)') '        <DataItem Dimensions="'//trim(itoa(nz))//' '//trim(itoa(ny))//' '//trim(itoa(nx))// &
           '" NumberType="Float" Precision="'//trim(itoa(int(bpr)))// &
           '" Format="Binary" Endian="Little" Seek="'//trim(i64toa(seek_bytes))// &
           '" >'//trim(fname)//'</DataItem>'
      write(u,'(a)') '      </Attribute>'
    end subroutine write_attr
    
    function itoa(i) result(s)
      implicit none
      integer, intent(in) :: i
      character(32) :: s
      write(s,'(i0)') i
    end function itoa

    function i64toa(i) result(s)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer(int64), intent(in) :: i
      character(64) :: s
      write(s,'(i0)') i
    end function i64toa

  end subroutine WriteXDMF
  
  pure function to_str(i) result(s)
    implicit none
    integer, intent(in) :: i
    character(len=32) :: s
    write(s,'(i0)') i
  end function to_str

end module mpiiomod
      
subroutine Output(is_final)
  use mpimod
  use basicmod
  use mpiiomod
  use config, only: asciiout
  implicit none

  logical, intent(in):: is_final

  real(8),save:: tout
  data tout / 0.0d0 /
  integer,save:: nout
  data nout / 1 /
  
  logical, save:: is_inited
  data is_inited /.false./
  
  logical,parameter:: binaryout= .true.


  if(time .lt. tout+dtout .and. .not. is_final) return

  call MPI_IO_Pack(nout)
  call MPI_IO_Write(nout)
  if(myid_w==0) call WriteXDMF(time,nout)
  
  if(asciiout) then
     call ASC_WRITE(nout)
  endif
  
  if(myid_w==0) print *, "output:",nout,time
          
  nout=nout+1
  tout=time
  
  return
  !         write(6,*) "bpf2",nflux2(mbps,i,j,k)
end subroutine Output


subroutine ASC_WRITE(nout)
  use mpimod
  use basicmod
  implicit none
  integer,intent(in):: nout
  integer:: i,j,k,n
  integer:: unitasc
  character(len=10),parameter :: dirname="ascdata/"
  character(40)::filename
  character(len=256):: header
  logical,save :: is_inited
  data is_inited / .false. /
  
  if (.not. is_inited) then
     call makedirs(dirname)
     is_inited =.true.
  endif
    
  write(filename,'(a4,i3.3,a1,i5.5,a4)')"snap",myid_w,"-",nout,".csv"
  filename = trim(dirname)//filename
    
  k=ks
  open(newunit=unitasc,file=filename,status='replace',form='formatted')
  write(unitasc,'((a1,1x),2(E12.3,1x))') "#",time,dt
  write(unitasc,'((a1,1x),i5)') "#",ie-is+1
  write(unitasc,'((a1,1x),i5)') "#",je-js+1
  write(unitasc,'((a1,1x),(i5,1x),E13.3)') "#",k,x3b(k)
  header = "# x y d vx vy p phi"
  do n=1,ncomp
     header = trim(header)//" X"//trim(adjustl(to_str(n)))
  enddo
  write(unitasc,'(A)') header 
  do j=js,je
     do i=is,ie
        write(unitasc,"(*(E13.3,1x))")x1b(i),x2b(j),d(i,j,k),v1(i,j,k),v2(i,j,k),p(i,j,k),gp(i,j,k),(Xcomp(n,i,j,k),n=1,ncomp)
     enddo
     write(unitasc,*) ""
  enddo
  close(unitasc)
    
  return
contains
  
  pure function to_str(i) result(s)
    implicit none
    integer, intent(in) :: i
    character(len=32) :: s
    write(s,'(i0)') i
  end function to_str

end subroutine ASC_WRITE




subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!  write(*, *) trim(command)
  call system(command)
end subroutine makedirs
