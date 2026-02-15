module boundarymod
  use config, only: periodicb, reflection, outflow &
       & , boundary_xin, boundary_xout &
       & , boundary_yin, boundary_yout &
       & , boundary_zin, boundary_zout &
       & , nomp
  use basicmod
  implicit none
  private
  
  integer,parameter:: nbc=9+ncomp
  !$omp declare target (nbc)
  integer,parameter:: nv1=3,nv2=4,nv3=5

  !-----------------------------------------------------------------
  ! GPU-offload friendly boundary buffers
  !-----------------------------------------------------------------

  logical :: bcbuf_ready = .false.

  real(8),allocatable :: varsendXstt(:,:,:,:), varsendXend(:,:,:,:)
  real(8),allocatable :: varsendYstt(:,:,:,:), varsendYend(:,:,:,:)
  real(8),allocatable :: varsendZstt(:,:,:,:), varsendZend(:,:,:,:)
  real(8),allocatable :: varrecvXstt(:,:,:,:), varrecvXend(:,:,:,:)
  real(8),allocatable :: varrecvYstt(:,:,:,:), varrecvYend(:,:,:,:)
  real(8),allocatable :: varrecvZstt(:,:,:,:), varrecvZend(:,:,:,:)

  !$omp declare target (varsendXstt,varsendXend,varsendYstt,varsendYend,varsendZstt,varsendZend)
  !$omp declare target (varrecvXstt,varrecvXend,varrecvYstt,varrecvYend,varrecvZstt,varrecvZend)
 
!!  real(8),dimension(mgn,jn,kn,nbc):: varsendXstt,varsendXend
!!  real(8),dimension(in,mgn,kn,nbc):: varsendYstt,varsendYend
!!  real(8),dimension(in,jn,mgn,nbc):: varsendZstt,varsendZend
!!  real(8),dimension(mgn,jn,kn,nbc):: varrecvXstt,varrecvXend
!!  real(8),dimension(in,mgn,kn,nbc):: varrecvYstt,varrecvYend
!!  real(8),dimension(in,jn,mgn,nbc):: varrecvZstt,varrecvZend
!!
!!!$omp declare target (varsendXstt,varsendXend)
!!!$omp declare target (varsendYstt,varsendYend)
!!!$omp declare target (varsendZstt,varsendZend)
!!!$omp declare target (varrecvXstt,varrecvXend)
!!!$omp declare target (varrecvYstt,varrecvYend)
!!!$omp declare target (varrecvZstt,varrecvZend)
!!  
  public:: BoundaryCondition
  public:: InitBoundaryBuffers, FinalizeBoundaryBuffers
contains

  subroutine InitBoundaryBuffers
    implicit none
    if (bcbuf_ready) return

    allocate(varsendXstt(mgn,jn,kn,nbc), varsendXend(mgn,jn,kn,nbc))
    allocate(varsendYstt(in ,mgn,kn,nbc), varsendYend(in ,mgn,kn,nbc))
    allocate(varsendZstt(in ,jn ,mgn,nbc), varsendZend(in ,jn ,mgn,nbc))
    allocate(varrecvXstt(mgn,jn,kn,nbc), varrecvXend(mgn,jn,kn,nbc))
    allocate(varrecvYstt(in ,mgn,kn,nbc), varrecvYend(in ,mgn,kn,nbc))
    allocate(varrecvZstt(in ,jn ,mgn,nbc), varrecvZend(in ,jn ,mgn,nbc))

    !$omp target enter data map(alloc: varsendXstt,varsendXend,varsendYstt,varsendYend,varsendZstt,varsendZend)
    !$omp target enter data map(alloc: varrecvXstt,varrecvXend,varrecvYstt,varrecvYend,varrecvZstt,varrecvZend)

    bcbuf_ready = .true.
  end subroutine InitBoundaryBuffers

  subroutine FinalizeBoundaryBuffers
    implicit none
    if (.not. bcbuf_ready) return

    !$omp target exit data map(delete: varsendXstt,varsendXend,varsendYstt,varsendYend,varsendZstt,varsendZend)
    !$omp target exit data map(delete: varrecvXstt,varrecvXend,varrecvYstt,varrecvYend,varrecvZstt,varrecvZend)

    deallocate(varsendXstt,varsendXend,varsendYstt,varsendYend,varsendZstt,varsendZend)
    deallocate(varrecvXstt,varrecvXend,varrecvYstt,varrecvYend,varrecvZstt,varrecvZend)

    bcbuf_ready = .false.
  end subroutine FinalizeBoundaryBuffers

  subroutine BoundaryCondition
    implicit none
    integer::i,j,k
    if (.not. bcbuf_ready) call InitBoundaryBuffers
  
!$omp target teams distribute parallel do collapse(3)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
     varsendXend(i,j,k,1) =  d(ie-mgn+i,j,k)
     varsendXend(i,j,k,2) = ei(ie-mgn+i,j,k)
     varsendXend(i,j,k,nv1) = v1(ie-mgn+i,j,k)
     varsendXend(i,j,k,nv2) = v2(ie-mgn+i,j,k)
     varsendXend(i,j,k,nv3) = v3(ie-mgn+i,j,k)
     varsendXend(i,j,k,6) = b1(ie-mgn+i,j,k)
     varsendXend(i,j,k,7) = b2(ie-mgn+i,j,k)
     varsendXend(i,j,k,8) = b3(ie-mgn+i,j,k)
     varsendXend(i,j,k,9) = bp(ie-mgn+i,j,k)
     varsendXend(i,j,k,10:nbc) = Xcomp(1:ncomp,ie-mgn+i,j,k)

     varsendXstt(i,j,k,1) =  d(  is+i-1,j,k)
     varsendXstt(i,j,k,2) = ei(  is+i-1,j,k)
     varsendXstt(i,j,k,nv1) = v1(  is+i-1,j,k)
     varsendXstt(i,j,k,nv2) = v2(  is+i-1,j,k)
     varsendXstt(i,j,k,nv3) = v3(  is+i-1,j,k)
     varsendXstt(i,j,k,6) = b1(  is+i-1,j,k)
     varsendXstt(i,j,k,7) = b2(  is+i-1,j,k)
     varsendXstt(i,j,k,8) = b3(  is+i-1,j,k)
     varsendXstt(i,j,k,9) = bp(  is+i-1,j,k)
     varsendXstt(i,j,k,10:nbc) = Xcomp(1:ncomp,is+i-1,j,k)
  enddo
  enddo
  enddo

!$omp target teams distribute parallel do collapse(3)
  do k=1,kn-1
  do i=1,in-1
  do j=1,mgn
     varsendYend(i,j,k,1) =  d(i,je-mgn+j,k)
     varsendYend(i,j,k,2) = ei(i,je-mgn+j,k)
     varsendYend(i,j,k,nv1) = v1(i,je-mgn+j,k)
     varsendYend(i,j,k,nv2) = v2(i,je-mgn+j,k)
     varsendYend(i,j,k,nv3) = v3(i,je-mgn+j,k)
     varsendYend(i,j,k,6) = b1(i,je-mgn+j,k)
     varsendYend(i,j,k,7) = b2(i,je-mgn+j,k)
     varsendYend(i,j,k,8) = b3(i,je-mgn+j,k)
     varsendYend(i,j,k,9) = bp(i,je-mgn+j,k)
     varsendYend(i,j,k,10:nbc) = Xcomp(1:ncomp,i,je-mgn+j,k)

     varsendYstt(i,j,k,1) =  d(i,  js+j-1,k)
     varsendYstt(i,j,k,2) = ei(i,  js+j-1,k)
     varsendYstt(i,j,k,nv1) = v1(i,  js+j-1,k)
     varsendYstt(i,j,k,nv2) = v2(i,  js+j-1,k)
     varsendYstt(i,j,k,nv3) = v3(i,  js+j-1,k)
     varsendYstt(i,j,k,6) = b1(i,  js+j-1,k)
     varsendYstt(i,j,k,7) = b2(i,  js+j-1,k)
     varsendYstt(i,j,k,8) = b3(i,  js+j-1,k)
     varsendYstt(i,j,k,9) = bp(i,  js+j-1,k)
     varsendYstt(i,j,k,10:nbc) = Xcomp(1:ncomp,i,  js+j-1,k)
  enddo
  enddo
  enddo

!$omp target teams distribute parallel do collapse(3)
  do j=1,jn-1
  !$omp target update from(varsendXstt,varsendXend,varsendYstt,varsendYend,varsendZstt,varsendZend)
  do i=1,in-1
  do k=1,mgn
     varsendZend(i,j,k,1) =  d(i,j,ke-mgn+k)
     varsendZend(i,j,k,2) = ei(i,j,ke-mgn+k)
     varsendZend(i,j,k,nv1) = v1(i,j,ke-mgn+k)
     varsendZend(i,j,k,nv2) = v2(i,j,ke-mgn+k)
     varsendZend(i,j,k,nv3) = v3(i,j,ke-mgn+k)
     varsendZend(i,j,k,6) = b1(i,j,ke-mgn+k)
     varsendZend(i,j,k,7) = b2(i,j,ke-mgn+k)
     varsendZend(i,j,k,8) = b3(i,j,ke-mgn+k)
     varsendZend(i,j,k,9) = bp(i,j,ke-mgn+k)
     varsendZend(i,j,k,10:nbc) = Xcomp(1:ncomp,i,j,ke-mgn+k)

     varsendZstt(i,j,k,1) =  d(i,j,ks+k-1  )
     varsendZstt(i,j,k,2) = ei(i,j,ks+k-1  )
     varsendZstt(i,j,k,nv1) = v1(i,j,ks+k-1  )
     varsendZstt(i,j,k,nv2) = v2(i,j,ks+k-1  )
     varsendZstt(i,j,k,nv3) = v3(i,j,ks+k-1  )
     varsendZstt(i,j,k,6) = b1(i,j,ks+k-1  )
     varsendZstt(i,j,k,7) = b2(i,j,ks+k-1  )
     varsendZstt(i,j,k,8) = b3(i,j,ks+k-1  )
     varsendZstt(i,j,k,9) = bp(i,j,ks+k-1  )
     varsendZstt(i,j,k,10:nbc) = Xcomp(1:ncomp,i,j,ks+k-1  )
  enddo
  enddo
  enddo

  call XbcSendRecv(varsendXstt,varsendXend,varrecvXstt,varrecvXend)
  call YbcSendRecv(varsendYstt,varsendYend,varrecvYstt,varrecvYend)
  call ZbcSendRecv(varsendZstt,varsendZend,varrecvZstt,varrecvZend)

  !$omp target update to(varrecvXstt,varrecvXend,varrecvYstt,varrecvYend,varrecvZstt,varrecvZend)
  
!$omp target teams distribute parallel do collapse(3)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
      d(i,j,k) = varrecvXstt(i,j,k,1)
     ei(i,j,k) = varrecvXstt(i,j,k,2)
     v1(i,j,k) = varrecvXstt(i,j,k,nv1)
     v2(i,j,k) = varrecvXstt(i,j,k,nv2)
     v3(i,j,k) = varrecvXstt(i,j,k,nv3)
     b1(i,j,k) = varrecvXstt(i,j,k,6)
     b2(i,j,k) = varrecvXstt(i,j,k,7)
     b3(i,j,k) = varrecvXstt(i,j,k,8)
     bp(i,j,k) = varrecvXstt(i,j,k,9)
     Xcomp(1:ncomp,i,j,k) = varrecvXstt(i,j,k,10:nbc)
     
      d(ie+i,j,k) = varrecvXend(i,j,k,1)
     ei(ie+i,j,k) = varrecvXend(i,j,k,2)
     v1(ie+i,j,k) = varrecvXend(i,j,k,nv1)
     v2(ie+i,j,k) = varrecvXend(i,j,k,nv2)
     v3(ie+i,j,k) = varrecvXend(i,j,k,nv3)
     b1(ie+i,j,k) = varrecvXend(i,j,k,6)
     b2(ie+i,j,k) = varrecvXend(i,j,k,7)
     b3(ie+i,j,k) = varrecvXend(i,j,k,8)
     bp(ie+i,j,k) = varrecvXend(i,j,k,9)
     Xcomp(1:ncomp,ie+i,j,k) = varrecvXend(i,j,k,10:nbc)
  enddo
  enddo
  enddo

!$omp target teams distribute parallel do collapse(3)
  do k=1,kn-1
  do i=1,in-1
  do j=1,mgn
      d(i,j,k) = varrecvYstt(i,j,k,1)
     ei(i,j,k) = varrecvYstt(i,j,k,2)
     v1(i,j,k) = varrecvYstt(i,j,k,nv1)
     v2(i,j,k) = varrecvYstt(i,j,k,nv2)
     v3(i,j,k) = varrecvYstt(i,j,k,nv3)
     b1(i,j,k) = varrecvYstt(i,j,k,6)
     b2(i,j,k) = varrecvYstt(i,j,k,7)
     b3(i,j,k) = varrecvYstt(i,j,k,8)
     bp(i,j,k) = varrecvYstt(i,j,k,9)
     Xcomp(1:ncomp,i,j,k) = varrecvYstt(i,j,k,10:nbc)
     
      d(i,je+j,k) = varrecvYend(i,j,k,1)
     ei(i,je+j,k) = varrecvYend(i,j,k,2)
     v1(i,je+j,k) = varrecvYend(i,j,k,nv1)
     v2(i,je+j,k) = varrecvYend(i,j,k,nv2)
     v3(i,je+j,k) = varrecvYend(i,j,k,nv3)
     b1(i,je+j,k) = varrecvYend(i,j,k,6)
     b2(i,je+j,k) = varrecvYend(i,j,k,7)
     b3(i,je+j,k) = varrecvYend(i,j,k,8)
     bp(i,je+j,k) = varrecvYend(i,j,k,9)
     Xcomp(1:ncomp,i,je+j,k) = varrecvYend(i,j,k,10:nbc)
  enddo
  enddo
  enddo


!$omp target teams distribute parallel do collapse(3)
  do j=1,jn-1
  do i=1,in-1
  do k=1,mgn
      d(i,j,k) = varrecvZstt(i,j,k,1)
     ei(i,j,k) = varrecvZstt(i,j,k,2)
     v1(i,j,k) = varrecvZstt(i,j,k,nv1)
     v2(i,j,k) = varrecvZstt(i,j,k,nv2)
     v3(i,j,k) = varrecvZstt(i,j,k,nv3)
     b1(i,j,k) = varrecvZstt(i,j,k,6)
     b2(i,j,k) = varrecvZstt(i,j,k,7)
     b3(i,j,k) = varrecvZstt(i,j,k,8)
     bp(i,j,k) = varrecvZstt(i,j,k,9)
     Xcomp(1:ncomp,i,j,k) = varrecvZstt(i,j,k,10:nbc)
     
      d(i,j,ke+k) = varrecvZend(i,j,k,1)
     ei(i,j,ke+k) = varrecvZend(i,j,k,2)
     v1(i,j,ke+k) = varrecvZend(i,j,k,nv1)
     v2(i,j,ke+k) = varrecvZend(i,j,k,nv2)
     v3(i,j,ke+k) = varrecvZend(i,j,k,nv3)
     b1(i,j,ke+k) = varrecvZend(i,j,k,6)
     b2(i,j,ke+k) = varrecvZend(i,j,k,7)
     b3(i,j,ke+k) = varrecvZend(i,j,k,8)
     bp(i,j,ke+k) = varrecvZend(i,j,k,9)
     Xcomp(1:ncomp,i,j,ke+k) = varrecvZend(i,j,k,10:nbc)
  enddo
  enddo
  enddo
  
  return
end subroutine BoundaryCondition

subroutine XbcSendRecv(varsendXstt,varsendXend,varrecvXstt,varrecvXend)
  use mpimod
  implicit none
  real(8),dimension(mgn,jn,kn,nbc),intent(in) ::varsendXstt,varsendXend
  real(8),dimension(mgn,jn,kn,nbc),intent(out)::varrecvXstt,varrecvXend
  integer::i,j,k,n
  
  single: if(ntiles(1) == 1) then
!$omp parallel do collapse(3)  private(n)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
     select case(boundary_xin)
     case(periodicb)
        do n=1,nbc
           varrecvXstt(i,j,k,n) = varsendXend(i,j,k,n)
        enddo
     case(reflection)
        do n=1,nbc
           varrecvXstt(i,j,k,n) = varsendXstt(mgn-i+1,j,k,n)
        enddo
        varrecvXstt(i,j,k,nv1) = -varrecvXstt(i,j,k,nv1)  ! v1
     case(outflow)
        do n=1,nbc
           varrecvXstt(i,j,k,n) = varsendXstt(1,j,k,n)
        enddo
     end select
     
     select case(boundary_xout)
     case(periodicb)
        do n=1,nbc
           varrecvXend(i,j,k,n) = varsendXstt(i,j,k,n)
        enddo
     case(reflection)
        do n=1,nbc
           varrecvXend(i,j,k,n) = varsendXend(mgn-i+1,j,k,n)
        enddo
        varrecvXend(i,j,k,nv1) = -varrecvXend(i,j,k,nv1)  ! v1
     case(outflow)
        do n=1,nbc
           varrecvXend(i,j,k,n) = varsendXend(mgn,j,k,n)
        enddo
     end select
     
  enddo
  enddo
  enddo
  
  else ! single

  n1mdir: if (n1m /= MPI_PROC_NULL) then
  else !n1dir

!$omp parallel do collapse(3)  private(n)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
     select case(boundary_xin)
     case(reflection)
        do n=1,nbc
           varrecvXstt(i,j,k,n) = varsendXstt(mgn-i+1,j,k,n)
        enddo
        varrecvXstt(i,j,k,nv1) = -varrecvXstt(i,j,k,nv1)  ! v1
     case(outflow)
        do n=1,nbc
           varrecvXstt(i,j,k,n) = varsendXstt(1,j,k,n)
        enddo
     end select
  enddo
  enddo
  enddo
     
  endif n1mdir
  
  n1pdir: if (n1p /= MPI_PROC_NULL) then
     
  else ! n1pdir
     
!$omp parallel do collapse(3)  private(n)
  do k=1,kn-1
  do j=1,jn-1
  do i=1,mgn
     select case(boundary_xout)
     case(reflection)
        do n=1,nbc
           varrecvXend(i,j,k,n) = varsendXend(mgn-i+1,j,k,n)
        enddo
        varrecvXend(i,j,k,nv1) = -varrecvXend(i,j,k,nv1)  ! v1
     case(outflow)
        do n=1,nbc
           varrecvXend(i,j,k,n) = varsendXend(mgn,j,k,n)
        enddo
     end select
  enddo
  enddo
  enddo
  endif n1pdir
  if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
  nreq = 0
  endif single
  
  
  return
end subroutine XbcSendRecv

subroutine YbcSendRecv(varsendYstt,varsendYend,varrecvYstt,varrecvYend)
  use mpimod
  implicit none
  real(8),dimension(in,mgn,kn,nbc),intent(in) ::varsendYstt,varsendYend
  real(8),dimension(in,mgn,kn,nbc),intent(out)::varrecvYstt,varrecvYend
  integer::i,j,k,n

  single: if(ntiles(2) == 1) then
!$omp parallel do collapse(3)  private(n)
  do k=1,kn-1
  do j=1,mgn
  do i=1,in-1
        
     select case(boundary_yin)
     case(periodicb)
        do n=1,nbc
           varrecvYstt(i,j,k,n) = varsendYend(i,j,k,n)
        enddo
     case(reflection)
        do n=1,nbc
           varrecvYstt(i,j,k,n) = varsendYstt(i,mgn-j+1,k,n)
        enddo
        varrecvYstt(i,j,k,nv2) = -varrecvYstt(i,j,k,nv2)  ! v2
     case(outflow)
        do n=1,nbc
           varrecvYstt(i,j,k,n) = varsendYstt(i,1,k,n)
        enddo
     end select
     
     select case(boundary_yout)
     case(periodicb)
        do n=1,nbc
           varrecvYend(i,j,k,n) = varsendYstt(i,j,k,n)
        enddo
     case(reflection)
        do n=1,nbc
           varrecvYend(i,j,k,n) = varsendYend(i,mgn-j+1,k,n)
        enddo
        varrecvYend(i,j,k,nv2) = -varrecvYend(i,j,k,nv2)  ! v2
     case(outflow)
        do n=1,nbc
           varrecvYend(i,j,k,n) = varsendYend(i,mgn,k,n)
        enddo
     end select

  enddo
  enddo
  enddo
  else ! single

     n2mdir: if (n2m /= MPI_PROC_NULL) then
     else
!$omp parallel do collapse(3)  private(n)
        do k=1,kn-1
        do j=1,mgn
        do i=1,in-1

           select case(boundary_yin)
           case(reflection)
              do n=1,nbc
                 varrecvYstt(i,j,k,n) = varsendYstt(i,mgn-j+1,k,n)
              enddo
              varrecvYstt(i,j,k,nv2) = -varrecvYstt(i,j,k,nv2)  ! v2
           case(outflow)
              do n=1,nbc
                 varrecvYstt(i,j,k,n) = varsendYstt(i,1,k,n)
              enddo
           end select
        enddo
        enddo
        enddo
     endif n2mdir
     
     n2pdir: if (n2p /= MPI_PROC_NULL) then
     else ! n2pdir
!$omp parallel do collapse(3)  private(n)
        do k=1,kn-1
        do i=1,in-1
        do j=1,mgn
           select case(boundary_yout)
           case(reflection)
              do n=1,nbc
                 varrecvYend(i,j,k,n) = varsendYend(i,mgn-j+1,k,n)
              enddo
              varrecvYend(i,j,k,nv2) = -varrecvYend(i,j,k,nv2)  ! v2
           case(outflow)
              do n=1,nbc
                 varrecvYend(i,j,k,n) = varsendYend(i,mgn,k,n)
              enddo
           end select
     
        enddo
        enddo
        enddo
        
        endif n2pdir
     
     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0
  endif single

  return
end subroutine YbcSendRecv

subroutine ZbcSendRecv(varsendZstt,varsendZend,varrecvZstt,varrecvZend)
  use mpimod
  implicit none
  real(8),dimension(in,jn,mgn,nbc),intent(in) ::varsendZstt,varsendZend
  real(8),dimension(in,jn,mgn,nbc),intent(out)::varrecvZstt,varrecvZend
  integer::i,j,k,n
  
  single: if(ntiles(3) == 1) then
!$omp parallel do collapse(3)  private(n)
  do k=1,mgn
  do j=1,jn-1
  do i=1,in-1
         
     select case(boundary_zin)
     case(periodicb)
        do n=1,nbc
           varrecvZstt(i,j,k,n) = varsendZend(i,j,k,n)
        enddo
     case(reflection)
        do n=1,nbc
           varrecvZstt(i,j,k,n) = varsendZstt(i,j,mgn-k+1,n)
        enddo
        varrecvZstt(i,j,k,nv3) = -varrecvZstt(i,j,k,nv3)  ! v3
     case(outflow)
        do n=1,nbc
           varrecvZstt(i,j,k,n) = varsendZstt(i,j,1,n)
        enddo
     end select
     
     select case(boundary_zout)
     case(periodicb)
        do n=1,nbc
           varrecvZend(i,j,k,n) = varsendZstt(i,j,k,n)
        enddo
     case(reflection)
        do n=1,nbc
           varrecvZend(i,j,k,n) = varsendZend(i,j,mgn-k+1,n)
        enddo
        varrecvZend(i,j,k,nv3) = -varrecvZend(i,j,k,nv3)  ! v3
     case(outflow)
        do n=1,nbc
           varrecvZend(i,j,k,n) = varsendZend(i,j,mgn,n)
        enddo
     end select

  enddo
  enddo
  enddo
  else ! single

     n3mdir: if (n3m /= MPI_PROC_NULL) then
     else
        
!$omp parallel do collapse(3)  private(n)
        do k=1,mgn
        do j=1,jn-1
        do i=1,in-1
         
           select case(boundary_zin)
           case(reflection)
              do n=1,nbc
                 varrecvZstt(i,j,k,n) = varsendZstt(i,j,mgn-k+1,n)
              enddo
              varrecvZstt(i,j,k,nv3) = -varrecvZstt(i,j,k,nv3)  ! v3
           case(outflow)
              do n=1,nbc
                 varrecvZstt(i,j,k,n) = varsendZstt(i,j,1,n)
              enddo
           end select
    
        enddo
        enddo
        enddo
     endif n3mdir
     
     n3pdir: if (n3p /= MPI_PROC_NULL) then
     else
!$omp parallel do collapse(3)  private(n)
        do j=1,jn-1
        do i=1,in-1
        do k=1,mgn
        
           select case(boundary_zout)
           case(reflection)
              do n=1,nbc
                 varrecvZend(i,j,k,n) = varsendZend(i,j,mgn-k+1,n)
              enddo
              varrecvZend(i,j,k,nv3) = -varrecvZend(i,j,k,nv3)  ! v3
           case(outflow)
              do n=1,nbc
                 varrecvZend(i,j,k,n) = varsendZend(i,j,mgn,n)
              enddo
           end select
           
  enddo
  enddo
  enddo
     endif n3pdir
     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0
  endif single

  return
end subroutine ZbcSendRecv
end  module boundarymod
