module mpimod
  use mpi
  implicit none
  integer, parameter :: mreq  = 300
  integer :: stat(MPI_STATUS_SIZE,mreq)                     
  integer :: req(mreq)
  
  integer :: ierr,myid_w, nprocs_w
  integer :: mpi_comm_hyd,myid_hyd, nprocs_hyd
  integer :: comm3d,myid, nprocs
  logical :: periodic(3)
  integer :: ntiles(3), coords(3)
  logical :: reorder
  integer :: n1m, n1p, n2m, n2p, n3m, n3p
  integer :: nreq, nsub
  integer ::   gpuid, ngpus
!$acc declare create(myid_w)
  
  real(8),dimension(2):: bufinpmin, bufoutmin
!$acc declare create(bufinpmin,bufoutmin)
  real(8),dimension(2):: bufinpmax, bufoutmax
!$acc declare create(bufinpmax,bufoutmax)

contains
subroutine InitializeMPI
  use openacc
  implicit none
  integer::key,color
  integer::np_hyd

! Initialize MPI
  call MPI_INIT( ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )
  
  ntiles(1)=1
  ntiles(2)=2
  ntiles(3)=2
  periodic(1)=.true.
  periodic(2)=.true.
  periodic(3)=.true.
  if(myid_w == 0) then
     print *, "MPI process=",nprocs_w
     print *, "decomposition=",ntiles(1),ntiles(2),ntiles(3)
  endif

  call MPI_BCAST(ntiles,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(periodic,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

! Making 3D strucure
  np_hyd = ntiles(1)*ntiles(2)*ntiles(3)
  color = int(myid_w/np_hyd)
  key   = myid_w   
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,mpi_comm_hyd,ierr)
  call MPI_COMM_SIZE( mpi_comm_hyd, nprocs_hyd, ierr )
  call MPI_COMM_RANK( mpi_comm_hyd, myid_hyd , ierr )     
  
! Create a virtual Cartesian topology for the domain decomposition.
!
  call MPI_CART_CREATE( mpi_comm_hyd, 3, ntiles, periodic &
       &                    , reorder, comm3d, ierr )
  call MPI_COMM_RANK( comm3d, myid,     ierr )
  call MPI_COMM_SIZE( comm3d, nprocs,   ierr )
!
! Find the ranks of my neighbors; find my virtual Cartesian coords.
!
  call MPI_CART_SHIFT( comm3d, 0, 1, n1m, n1p, ierr )
  call MPI_CART_SHIFT( comm3d, 1, 1, n2m, n2p, ierr )
  call MPI_CART_SHIFT( comm3d, 2, 1, n3m, n3p, ierr )
  !
  call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )

!> debug  
  call MPI_Comm_set_errhandler(comm3d, MPI_ERRORS_RETURN, ierr)
  
  ngpus = acc_get_num_devices(acc_device_nvidia)
  if(myid_w == 0) then
     print *, "num of GPUs = ", ngpus
  end if

  if(ngpus == 0) then
     gpuid = -1
  else
     gpuid = mod(myid_w, ngpus)
  endif
  
  if(gpuid >= 0) then
     call acc_set_device_num(gpuid, acc_device_nvidia)
  end if
!$acc update device (myid_w)
  return
end subroutine InitializeMPI

subroutine FinalizeMPI
  implicit none
  call MPI_FINALIZE(ierr)
end subroutine FinalizeMPI

subroutine MPIminfind
  implicit none
  integer :: err_len
  character(len=MPI_MAX_ERROR_STRING) :: err_string
!$acc host_data use_device(bufinpmin,bufoutmin)
       call MPI_ALLREDUCE( bufinpmin(1), bufoutmin(1), 1 &
     &                   , MPI_2DOUBLE_PRECISION   &
     &                   , MPI_MINLOC, comm3d, ierr)      
!$acc end host_data
       if (ierr /= MPI_SUCCESS) then
          call MPI_Error_string(ierr, err_string, err_len, ierr)
          print *,"error in MPIminfind", trim(err_string)
       endif
end subroutine MPIminfind

subroutine MPImaxfind
  implicit none
  integer :: err_len
  character(len=MPI_MAX_ERROR_STRING) :: err_string
!$acc host_data use_device(bufinpmax,bufoutmax)
       call MPI_ALLREDUCE( bufinpmax(1), bufoutmax(1), 1 &
     &                   , MPI_2DOUBLE_PRECISION   &
     &                   , MPI_MAXLOC, comm3d, ierr)
!$acc end host_data
       if (ierr /= MPI_SUCCESS) then
          call MPI_Error_string(ierr, err_string, err_len, ierr)
          print *,"error in MPIminfind", trim(err_string)
       endif
end subroutine MPImaxfind

end module mpimod

