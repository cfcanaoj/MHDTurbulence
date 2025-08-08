module boundarymod
  use mpimod
  use basicmod
  implicit none
  private
  real(8),dimension(mgn,jn,kn,nbc):: varsendXstt,varsendXend
  real(8),dimension(in,mgn,kn,nbc):: varsendYstt,varsendYend
  real(8),dimension(in,jn,mgn,nbc):: varsendZstt,varsendZend
  real(8),dimension(mgn,jn,kn,nbc):: varrecvXstt,varrecvXend
  real(8),dimension(in,mgn,kn,nbc):: varrecvYstt,varrecvYend
  real(8),dimension(in,jn,mgn,nbc):: varrecvZstt,varrecvZend

!$omp declare target(varsendXstt,varsendXend)
!$omp declare target(varsendYstt,varsendYend)
!$omp declare target(varsendZstt,varsendZend)
!$omp declare target(varrecvXstt,varrecvXend)
!$omp declare target(varrecvYstt,varrecvYend)
!$omp declare target(varrecvZstt,varrecvZend)
  
  public:: BoundaryCondition
contains
  subroutine BoundaryCondition
    implicit none
    integer::i,j,k
!$omp target defaultmap(none)
!$omp loop order(concurrent) collapse(3)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
     varsendXend(i,j,k,1) =  d(ie-mgn+i,j,k)
     varsendXend(i,j,k,2) = ei(ie-mgn+i,j,k)
     varsendXend(i,j,k,3) = v1(ie-mgn+i,j,k)
     varsendXend(i,j,k,4) = v2(ie-mgn+i,j,k)
     varsendXend(i,j,k,5) = v3(ie-mgn+i,j,k)
     varsendXend(i,j,k,6) = b1(ie-mgn+i,j,k)
     varsendXend(i,j,k,7) = b2(ie-mgn+i,j,k)
     varsendXend(i,j,k,8) = b3(ie-mgn+i,j,k)
     varsendXend(i,j,k,9) = bp(ie-mgn+i,j,k)

     varsendXstt(i,j,k,1) =  d(  is+i-1,j,k)
     varsendXstt(i,j,k,2) = ei(  is+i-1,j,k)
     varsendXstt(i,j,k,3) = v1(  is+i-1,j,k)
     varsendXstt(i,j,k,4) = v2(  is+i-1,j,k)
     varsendXstt(i,j,k,5) = v3(  is+i-1,j,k)
     varsendXstt(i,j,k,6) = b1(  is+i-1,j,k)
     varsendXstt(i,j,k,7) = b2(  is+i-1,j,k)
     varsendXstt(i,j,k,8) = b3(  is+i-1,j,k)
     varsendXstt(i,j,k,9) = bp(  is+i-1,j,k)
  enddo
  enddo
  enddo

!$omp loop order(concurrent) collapse(3)
  do k=1,kn-1
  do i=1,in-1
  do j=1,mgn
     varsendYend(i,j,k,1) =  d(i,je-mgn+j,k)
     varsendYend(i,j,k,2) = ei(i,je-mgn+j,k)
     varsendYend(i,j,k,3) = v1(i,je-mgn+j,k)
     varsendYend(i,j,k,4) = v2(i,je-mgn+j,k)
     varsendYend(i,j,k,5) = v3(i,je-mgn+j,k)
     varsendYend(i,j,k,6) = b1(i,je-mgn+j,k)
     varsendYend(i,j,k,7) = b2(i,je-mgn+j,k)
     varsendYend(i,j,k,8) = b3(i,je-mgn+j,k)
     varsendYend(i,j,k,9) = bp(i,je-mgn+j,k)

     varsendYstt(i,j,k,1) =  d(i,  js+j-1,k)
     varsendYstt(i,j,k,2) = ei(i,  js+j-1,k)
     varsendYstt(i,j,k,3) = v1(i,  js+j-1,k)
     varsendYstt(i,j,k,4) = v2(i,  js+j-1,k)
     varsendYstt(i,j,k,5) = v3(i,  js+j-1,k)
     varsendYstt(i,j,k,6) = b1(i,  js+j-1,k)
     varsendYstt(i,j,k,7) = b2(i,  js+j-1,k)
     varsendYstt(i,j,k,8) = b3(i,  js+j-1,k)
     varsendYstt(i,j,k,9) = bp(i,  js+j-1,k)
  enddo
  enddo
  enddo

!$omp loop order(concurrent) collapse(3)
  do j=1,jn-1
  do i=1,in-1
  do k=1,mgn
     varsendZend(i,j,k,1) =  d(i,j,ke-mgn+k)
     varsendZend(i,j,k,2) = ei(i,j,ke-mgn+k)
     varsendZend(i,j,k,3) = v1(i,j,ke-mgn+k)
     varsendZend(i,j,k,4) = v2(i,j,ke-mgn+k)
     varsendZend(i,j,k,5) = v3(i,j,ke-mgn+k)
     varsendZend(i,j,k,6) = b1(i,j,ke-mgn+k)
     varsendZend(i,j,k,8) = b3(i,j,ke-mgn+k)
     varsendZend(i,j,k,9) = bp(i,j,ke-mgn+k)

     varsendZstt(i,j,k,1) =  d(i,j,ks+k-1  )
     varsendZstt(i,j,k,2) = ei(i,j,ks+k-1  )
     varsendZstt(i,j,k,3) = v1(i,j,ks+k-1  )
     varsendZstt(i,j,k,4) = v2(i,j,ks+k-1  )
     varsendZstt(i,j,k,5) = v3(i,j,ks+k-1  )
     varsendZstt(i,j,k,6) = b1(i,j,ks+k-1  )
     varsendZstt(i,j,k,7) = b2(i,j,ks+k-1  )
     varsendZstt(i,j,k,8) = b3(i,j,ks+k-1  )
     varsendZstt(i,j,k,9) = bp(i,j,ks+k-1  )
  enddo
  enddo
  enddo
!$omp end target


  call XbcSendRecv
  call YbcSendRecv
  call ZbcSendRecv
  
!$omp target defaultmap(none)
!$omp loop order(concurrent) collapse(3)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
      d(i,j,k) = varrecvXstt(i,j,k,1)
     ei(i,j,k) = varrecvXstt(i,j,k,2)
     v1(i,j,k) = varrecvXstt(i,j,k,3)
     v2(i,j,k) = varrecvXstt(i,j,k,4)
     v3(i,j,k) = varrecvXstt(i,j,k,5)
     b1(i,j,k) = varrecvXstt(i,j,k,6)
     b2(i,j,k) = varrecvXstt(i,j,k,7)
     b3(i,j,k) = varrecvXstt(i,j,k,8)
     bp(i,j,k) = varrecvXstt(i,j,k,9)
     
      d(ie+i,j,k) = varrecvXend(i,j,k,1)
     ei(ie+i,j,k) = varrecvXend(i,j,k,2)
     v1(ie+i,j,k) = varrecvXend(i,j,k,3)
     v2(ie+i,j,k) = varrecvXend(i,j,k,4)
     v3(ie+i,j,k) = varrecvXend(i,j,k,5)
     b1(ie+i,j,k) = varrecvXend(i,j,k,6)
     b2(ie+i,j,k) = varrecvXend(i,j,k,7)
     b3(ie+i,j,k) = varrecvXend(i,j,k,8)
     bp(ie+i,j,k) = varrecvXend(i,j,k,9)
  enddo
  enddo
  enddo

!$omp loop order(concurrent) collapse(3)
  do k=1,kn-1
  do i=1,in-1
  do j=1,mgn
      d(i,j,k) = varrecvYstt(i,j,k,1)
     ei(i,j,k) = varrecvYstt(i,j,k,2)
     v1(i,j,k) = varrecvYstt(i,j,k,3)
     v2(i,j,k) = varrecvYstt(i,j,k,4)
     v3(i,j,k) = varrecvYstt(i,j,k,5)
     b1(i,j,k) = varrecvYstt(i,j,k,6)
     b2(i,j,k) = varrecvYstt(i,j,k,7)
     b3(i,j,k) = varrecvYstt(i,j,k,8)
     bp(i,j,k) = varrecvYstt(i,j,k,9)
     
      d(i,je+j,k) = varrecvYend(i,j,k,1)
     ei(i,je+j,k) = varrecvYend(i,j,k,2)
     v1(i,je+j,k) = varrecvYend(i,j,k,3)
     v2(i,je+j,k) = varrecvYend(i,j,k,4)
     v3(i,je+j,k) = varrecvYend(i,j,k,5)
     b1(i,je+j,k) = varrecvYend(i,j,k,6)
     b2(i,je+j,k) = varrecvYend(i,j,k,7)
     b3(i,je+j,k) = varrecvYend(i,j,k,8)
     bp(i,je+j,k) = varrecvYend(i,j,k,9)
  enddo
  enddo
  enddo


!$omp loop order(concurrent) collapse(3)
  do j=1,jn-1
  do i=1,in-1
  do k=1,mgn
      d(i,j,k) = varrecvZstt(i,j,k,1)
     ei(i,j,k) = varrecvZstt(i,j,k,2)
     v1(i,j,k) = varrecvZstt(i,j,k,3)
     v2(i,j,k) = varrecvZstt(i,j,k,4)
     v3(i,j,k) = varrecvZstt(i,j,k,5)
     b1(i,j,k) = varrecvZstt(i,j,k,6)
     b2(i,j,k) = varrecvZstt(i,j,k,7)
     b3(i,j,k) = varrecvZstt(i,j,k,8)
     bp(i,j,k) = varrecvZstt(i,j,k,9)
     
      d(i,j,ke+k) = varrecvZend(i,j,k,1)
     ei(i,j,ke+k) = varrecvZend(i,j,k,2)
     v1(i,j,ke+k) = varrecvZend(i,j,k,3)
     v2(i,j,ke+k) = varrecvZend(i,j,k,4)
     v3(i,j,ke+k) = varrecvZend(i,j,k,5)
     b1(i,j,ke+k) = varrecvZend(i,j,k,6)
     b2(i,j,ke+k) = varrecvZend(i,j,k,7)
     b3(i,j,ke+k) = varrecvZend(i,j,k,8)
     bp(i,j,ke+k) = varrecvZend(i,j,k,9)
  enddo
  enddo
  enddo
!$omp end target
  
  return
end subroutine BoundaryCondition

subroutine XbcSendRecv
  implicit none
  integer::i,j,k,n
  if(ntiles(1) == 1) then
!$omp target defaultmap(none)
!$omp loop order(concurrent) collapse(4)
  do n=1,nbc
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
     varrecvXstt(i,j,k,n) = varsendXend(i,j,k,n)
     varrecvXend(i,j,k,n) = varsendXstt(i,j,k,n)
  enddo
  enddo
  enddo
  enddo
!$omp end target
  
  else

!$omp target update from(varsendXstt,varsendXend)
     nreq = nreq + 1         
     call MPI_IRECV(varrecvXstt,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1m,1100, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendXstt,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1m, 1200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_IRECV(varrecvXend,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1p,1200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendXend,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1p, 1100, comm3d, req(nreq), ierr)

     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0
!$omp target update to(varrecvXstt,varrecvXend)

  endif

  return
end subroutine XbcSendRecv

subroutine YbcSendRecv
  implicit none
  integer::i,j,k,n
  if(ntiles(2) == 1) then
!$omp target defaultmap(none)
!$omp loop order(concurrent) collapse(4)
  do n=1,nbc
  do k=1,kn-1
  do j=1,mgn
  do i=1,in-1
     varrecvYstt(i,j,k,n) = varsendYend(i,j,k,n)
     varrecvYend(i,j,k,n) = varsendYstt(i,j,k,n)
  enddo
  enddo
  enddo
  enddo
!$omp end target
  else

!$omp target update from(varsendYstt,varsendYend)
     nreq = nreq + 1         
     call MPI_IRECV(varrecvYstt,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2m, 2100, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendYstt,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2m, 2200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_IRECV(varrecvYend,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2p,2200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendYend,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2p, 2100, comm3d, req(nreq), ierr)

     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0
!$omp target update to(varrecvYstt,varrecvYend)
  endif

  return
end subroutine YbcSendRecv

subroutine ZbcSendRecv
  implicit none
  integer::i,j,k,n
  
  if(ntiles(3) == 1) then
!$omp target defaultmap(none)
!$omp loop order(concurrent) collapse(4)
  do n=1,nbc
  do k=1,mgn
  do j=1,jn-1
  do i=1,in-1
     varrecvZstt(i,j,k,n) = varsendZend(i,j,k,n)
     varrecvZend(i,j,k,n) = varsendZstt(i,j,k,n)
  enddo
  enddo
  enddo
  enddo
!$omp end target
  else

!$omp target update from(varsendZstt,varsendZend)
     nreq = nreq + 1         
     call MPI_IRECV(varrecvZstt,mgn*in*jn*nbc &
    & , MPI_DOUBLE &
    & , n3m, 3100, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendZstt,mgn*in*jn*nbc &
    & , MPI_DOUBLE &
    & , n3m, 3200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_IRECV(varrecvZend,mgn*in*jn*nbc &
    & , MPI_DOUBLE &
    & , n3p, 3200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendZend,mgn*in*jn*nbc &
    & , MPI_DOUBLE &
    & , n3p, 3100, comm3d, req(nreq), ierr)

     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0
!$omp target update to(varrecvZstt,varrecvZend)
  endif

  return
end subroutine ZbcSendRecv
end  module boundarymod

! Code was translated using: /gwork0/takiwkkz/MHDturbArxiv/acctoomp/src/intel-application-migration-tool-for-openacc-to-openmp -suppress-openacc boundary.f90
