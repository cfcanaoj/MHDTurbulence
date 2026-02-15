module boundarymod
  use config, only: periodicb, reflection, outflow, &
       & boundary_xin, boundary_xout, &
       & boundary_yin, boundary_yout, &
       & boundary_zin, boundary_zout
  use basicmod
  implicit none

  integer,parameter :: nv1 = 3, nv2 = 4, nv3 = 5
  integer,parameter :: nbc = 9 + ncomp
  
  real(8),dimension(mgn,jn,kn,nbc):: BsXstt,BsXend
  real(8),dimension(in,mgn,kn,nbc):: BsYstt,BsYend
  real(8),dimension(in,jn,mgn,nbc):: BsZstt,BsZend
  real(8),dimension(mgn,jn,kn,nbc):: BrXstt,BrXend
  real(8),dimension(in,mgn,kn,nbc):: BrYstt,BrYend
  real(8),dimension(in,jn,mgn,nbc):: BrZstt,BrZend
    
!$omp declare target (BsXtt,BsXend)
!$omp declare target (BsYtt,BsYend)
!$omp declare target (BsZtt,BsZend)
!$omp declare target (BrXtt,BrXend)
!$omp declare target (BrYtt,BrYend)
!$omp declare target (BrZtt,BrZend)

contains

  subroutine BoundaryCondition
    implicit none
    integer :: i,j,k

!$omp target teams distribute parallel do collapse(3)
    do k=1,kn-1
    do j=1,jn-1
    do i=1,mgn
      BsXend(i,j,k,1) =  d(ie-mgn+i,j,k)
      BsXend(i,j,k,2) = ei(ie-mgn+i,j,k)
      BsXend(i,j,k,nv1) = v1(ie-mgn+i,j,k)
      BsXend(i,j,k,nv2) = v2(ie-mgn+i,j,k)
      BsXend(i,j,k,nv3) = v3(ie-mgn+i,j,k)
      BsXend(i,j,k,6) = b1(ie-mgn+i,j,k)
      BsXend(i,j,k,7) = b2(ie-mgn+i,j,k)
      BsXend(i,j,k,8) = b3(ie-mgn+i,j,k)
      BsXend(i,j,k,9) = bp(ie-mgn+i,j,k)
      BsXend(i,j,k,10:nbc) = Xcomp(1:ncomp,ie-mgn+i,j,k)

      BsXstt(i,j,k,1) =  d(is+i-1,j,k)
      BsXstt(i,j,k,2) = ei(is+i-1,j,k)
      BsXstt(i,j,k,nv1) = v1(is+i-1,j,k)
      BsXstt(i,j,k,nv2) = v2(is+i-1,j,k)
      BsXstt(i,j,k,nv3) = v3(is+i-1,j,k)
      BsXstt(i,j,k,6) = b1(is+i-1,j,k)
      BsXstt(i,j,k,7) = b2(is+i-1,j,k)
      BsXstt(i,j,k,8) = b3(is+i-1,j,k)
      BsXstt(i,j,k,9) = bp(is+i-1,j,k)
      BsXstt(i,j,k,10:nbc) = Xcomp(1:ncomp,is+i-1,j,k)
    end do
    end do
    end do

!$omp target teams distribute parallel do collapse(3)
    do k=1,kn-1
    do i=1,in-1
    do j=1,mgn
      BsYend(i,j,k,1) =  d(i,je-mgn+j,k)
      BsYend(i,j,k,2) = ei(i,je-mgn+j,k)
      BsYend(i,j,k,nv1) = v1(i,je-mgn+j,k)
      BsYend(i,j,k,nv2) = v2(i,je-mgn+j,k)
      BsYend(i,j,k,nv3) = v3(i,je-mgn+j,k)
      BsYend(i,j,k,6) = b1(i,je-mgn+j,k)
      BsYend(i,j,k,7) = b2(i,je-mgn+j,k)
      BsYend(i,j,k,8) = b3(i,je-mgn+j,k)
      BsYend(i,j,k,9) = bp(i,je-mgn+j,k)
      BsYend(i,j,k,10:nbc) = Xcomp(1:ncomp,i,je-mgn+j,k)

      BsYstt(i,j,k,1) =  d(i,js+j-1,k)
      BsYstt(i,j,k,2) = ei(i,js+j-1,k)
      BsYstt(i,j,k,nv1) = v1(i,js+j-1,k)
      BsYstt(i,j,k,nv2) = v2(i,js+j-1,k)
      BsYstt(i,j,k,nv3) = v3(i,js+j-1,k)
      BsYstt(i,j,k,6) = b1(i,js+j-1,k)
      BsYstt(i,j,k,7) = b2(i,js+j-1,k)
      BsYstt(i,j,k,8) = b3(i,js+j-1,k)
      BsYstt(i,j,k,9) = bp(i,js+j-1,k)
      BsYstt(i,j,k,10:nbc) = Xcomp(1:ncomp,i,js+j-1,k)
    end do
    end do
    end do

!$omp target teams distribute parallel do collapse(3)
    do j=1,jn-1
    do i=1,in-1
    do k=1,mgn
      BsZend(i,j,k,1) =  d(i,j,ke-mgn+k)
      BsZend(i,j,k,2) = ei(i,j,ke-mgn+k)
      BsZend(i,j,k,nv1) = v1(i,j,ke-mgn+k)
      BsZend(i,j,k,nv2) = v2(i,j,ke-mgn+k)
      BsZend(i,j,k,nv3) = v3(i,j,ke-mgn+k)
      BsZend(i,j,k,6) = b1(i,j,ke-mgn+k)
      BsZend(i,j,k,7) = b2(i,j,ke-mgn+k)
      BsZend(i,j,k,8) = b3(i,j,ke-mgn+k)
      BsZend(i,j,k,9) = bp(i,j,ke-mgn+k)
      BsZend(i,j,k,10:nbc) = Xcomp(1:ncomp,i,j,ke-mgn+k)

      BsZstt(i,j,k,1) =  d(i,j,ks+k-1)
      BsZstt(i,j,k,2) = ei(i,j,ks+k-1)
      BsZstt(i,j,k,nv1) = v1(i,j,ks+k-1)
      BsZstt(i,j,k,nv2) = v2(i,j,ks+k-1)
      BsZstt(i,j,k,nv3) = v3(i,j,ks+k-1)
      BsZstt(i,j,k,6) = b1(i,j,ks+k-1)
      BsZstt(i,j,k,7) = b2(i,j,ks+k-1)
      BsZstt(i,j,k,8) = b3(i,j,ks+k-1)
      BsZstt(i,j,k,9) = bp(i,j,ks+k-1)
      BsZstt(i,j,k,10:nbc) = Xcomp(1:ncomp,i,j,ks+k-1)
    end do
    end do
    end do

    ! ----------------------
    ! Exchange / apply physical BCs (MPI on host, BC possibly on device)
    ! ----------------------
    call XbcSendRecv()
    call YbcSendRecv()
    call ZbcSendRecv()

    ! ----------------------
    ! Unpack from receive buffers on device
    ! ----------------------

!$omp target teams distribute parallel do collapse(3)
    do k=1,kn-1
    do j=1,jn-1
    do i=1,mgn
       d(i,j,k) = BrXstt(i,j,k,1)
      ei(i,j,k) = BrXstt(i,j,k,2)
      v1(i,j,k) = BrXstt(i,j,k,nv1)
      v2(i,j,k) = BrXstt(i,j,k,nv2)
      v3(i,j,k) = BrXstt(i,j,k,nv3)
      b1(i,j,k) = BrXstt(i,j,k,6)
      b2(i,j,k) = BrXstt(i,j,k,7)
      b3(i,j,k) = BrXstt(i,j,k,8)
      bp(i,j,k) = BrXstt(i,j,k,9)
      Xcomp(1:ncomp,i,j,k) = BrXstt(i,j,k,10:nbc)

       d(ie+i,j,k) = BrXend(i,j,k,1)
      ei(ie+i,j,k) = BrXend(i,j,k,2)
      v1(ie+i,j,k) = BrXend(i,j,k,nv1)
      v2(ie+i,j,k) = BrXend(i,j,k,nv2)
      v3(ie+i,j,k) = BrXend(i,j,k,nv3)
      b1(ie+i,j,k) = BrXend(i,j,k,6)
      b2(ie+i,j,k) = BrXend(i,j,k,7)
      b3(ie+i,j,k) = BrXend(i,j,k,8)
      bp(ie+i,j,k) = BrXend(i,j,k,9)
      Xcomp(1:ncomp,ie+i,j,k) = BrXend(i,j,k,10:nbc)
    end do
    end do
    end do

!$omp target teams distribute parallel do collapse(3)
    do k=1,kn-1
    do j=1,mgn
    do i=1,in-1
       d(i,j,k) = BrYstt(i,j,k,1)
      ei(i,j,k) = BrYstt(i,j,k,2)
      v1(i,j,k) = BrYstt(i,j,k,nv1)
      v2(i,j,k) = BrYstt(i,j,k,nv2)
      v3(i,j,k) = BrYstt(i,j,k,nv3)
      b1(i,j,k) = BrYstt(i,j,k,6)
      b2(i,j,k) = BrYstt(i,j,k,7)
      b3(i,j,k) = BrYstt(i,j,k,8)
      bp(i,j,k) = BrYstt(i,j,k,9)
      Xcomp(1:ncomp,i,j,k) = BrYstt(i,j,k,10:nbc)

       d(i,je+j,k) = BrYend(i,j,k,1)
      ei(i,je+j,k) = BrYend(i,j,k,2)
      v1(i,je+j,k) = BrYend(i,j,k,nv1)
      v2(i,je+j,k) = BrYend(i,j,k,nv2)
      v3(i,je+j,k) = BrYend(i,j,k,nv3)
      b1(i,je+j,k) = BrYend(i,j,k,6)
      b2(i,je+j,k) = BrYend(i,j,k,7)
      b3(i,je+j,k) = BrYend(i,j,k,8)
      bp(i,je+j,k) = BrYend(i,j,k,9)
      Xcomp(1:ncomp,i,je+j,k) = BrYend(i,j,k,10:nbc)
    end do
    end do
    end do

!$omp target teams distribute parallel do collapse(3)
    do k=1,mgn
    do j=1,jn-1
    do i=1,in-1
       d(i,j,k) = BrZstt(i,j,k,1)
      ei(i,j,k) = BrZstt(i,j,k,2)
      v1(i,j,k) = BrZstt(i,j,k,nv1)
      v2(i,j,k) = BrZstt(i,j,k,nv2)
      v3(i,j,k) = BrZstt(i,j,k,nv3)
      b1(i,j,k) = BrZstt(i,j,k,6)
      b2(i,j,k) = BrZstt(i,j,k,7)
      b3(i,j,k) = BrZstt(i,j,k,8)
      bp(i,j,k) = BrZstt(i,j,k,9)
      Xcomp(1:ncomp,i,j,k) = BrZstt(i,j,k,10:nbc)

       d(i,j,ke+k) = BrZend(i,j,k,1)
      ei(i,j,ke+k) = BrZend(i,j,k,2)
      v1(i,j,ke+k) = BrZend(i,j,k,nv1)
      v2(i,j,ke+k) = BrZend(i,j,k,nv2)
      v3(i,j,ke+k) = BrZend(i,j,k,nv3)
      b1(i,j,ke+k) = BrZend(i,j,k,6)
      b2(i,j,ke+k) = BrZend(i,j,k,7)
      b3(i,j,ke+k) = BrZend(i,j,k,8)
      bp(i,j,ke+k) = BrZend(i,j,k,9)
      Xcomp(1:ncomp,i,j,ke+k) = BrZend(i,j,k,10:nbc)
    end do
    end do
    end do

  end subroutine BoundaryCondition


  subroutine XbcSendRecv()
    use mpimod
    implicit none
    integer :: i,j,k,n

    single: if (ntiles(1) == 1) then
      ! Local-only: compute halos on device
!$omp target teams distribute parallel do collapse(3)
      do k=1,kn-1
      do j=1,jn-1
      do i=1,mgn
        select case(boundary_xin)
        case(periodicb)
          BrXstt(i,j,k,1:nbc) = BsXend(i,j,k,1:nbc)
        case(reflection)
          BrXstt(i,j,k,1:nbc) = BsXstt(mgn-i+1,j,k,1:nbc)
          BrXstt(i,j,k,nv1) = -BrXstt(i,j,k,nv1)
        case(outflow)
          BrXstt(i,j,k,1:nbc) = BsXstt(1,j,k,1:nbc)
        end select

        select case(boundary_xout)
        case(periodicb)
          BrXend(i,j,k,1:nbc) = BsXstt(i,j,k,1:nbc)
        case(reflection)
          BrXend(i,j,k,1:nbc) = BsXend(mgn-i+1,j,k,1:nbc)
          BrXend(i,j,k,nv1) = -BrXend(i,j,k,nv1)
        case(outflow)
          BrXend(i,j,k,1:nbc) = BsXend(mgn,j,k,1:nbc)
        end select
      end do
      end do
      end do
   
      return
    end if single
    
    ! Ensure host has latest send buffers
!$omp target update from(BsXstt, BsXend)

    ! Post receives/sends where neighbors exist. Physical BC when MPI_PROC_NULL.
    if (n1m /= MPI_PROC_NULL) then
      call MPI_IRECV(BrXstt, size(BrXstt), MPI_DOUBLE_PRECISION, n1m, 1100, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
      call MPI_ISEND(BsXstt, size(BsXstt), MPI_DOUBLE_PRECISION, n1m, 1200, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
    else
      if (boundary_xin == reflection) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,jn-1; do i=1,mgn
          BrXstt(i,j,k,1:nbc) = BsXstt(mgn-i+1,j,k,1:nbc)
          BrXstt(i,j,k,nv1) = -BrXstt(i,j,k,nv1)
        end do; end do; end do
      else if (boundary_xin == outflow) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,jn-1; do i=1,mgn
          BrXstt(i,j,k,1:nbc) = BsXstt(1,j,k,1:nbc)
        end do; end do; end do
      end if
    end if

    if (n1p /= MPI_PROC_NULL) then
      call MPI_IRECV(BrXend, size(BrXend), MPI_DOUBLE_PRECISION, n1p, 1200, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
      call MPI_ISEND(BsXend, size(BsXend), MPI_DOUBLE_PRECISION, n1p, 1100, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
    else
      if (boundary_xout == reflection) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,jn-1; do i=1,mgn
          BrXend(i,j,k,1:nbc) = BsXend(mgn-i+1,j,k,1:nbc)
          BrXend(i,j,k,nv1) = -BrXend(i,j,k,nv1)
        end do; end do; end do
      else if (boundary_xout == outflow) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,jn-1; do i=1,mgn
          BrXend(i,j,k,1:nbc) = BsXend(mgn,j,k,1:nbc)
        end do; end do; end do
      end if
    end if

    if (nreq /= 0) then
      call MPI_WAITALL(nreq, req, stat, ierr)
      ! Push only faces that came from MPI
      if (n1m /= MPI_PROC_NULL) then
!$omp target update to(BrXstt)
      end if
      if (n1p /= MPI_PROC_NULL) then
!$omp target update to(BrXend)
      end if
      nreq = 0
    end if

  end subroutine XbcSendRecv


  subroutine YbcSendRecv()
    use mpimod
    implicit none
    integer :: i,j,k

    if (ntiles(2) == 1) then
!$omp target teams distribute parallel do collapse(3)
      do k=1,kn-1
      do j=1,mgn
      do i=1,in-1
        select case(boundary_yin)
        case(periodicb)
          BrYstt(i,j,k,1:nbc) = BsYend(i,j,k,1:nbc)
        case(reflection)
          BrYstt(i,j,k,1:nbc) = BsYstt(i,mgn-j+1,k,1:nbc)
          BrYstt(i,j,k,nv2) = -BrYstt(i,j,k,nv2)
        case(outflow)
          BrYstt(i,j,k,1:nbc) = BsYstt(i,1,k,1:nbc)
        end select

        select case(boundary_yout)
        case(periodicb)
          BrYend(i,j,k,1:nbc) = BsYstt(i,j,k,1:nbc)
        case(reflection)
          BrYend(i,j,k,1:nbc) = BsYend(i,mgn-j+1,k,1:nbc)
          BrYend(i,j,k,nv2) = -BrYend(i,j,k,nv2)
        case(outflow)
          BrYend(i,j,k,1:nbc) = BsYend(i,mgn,k,1:nbc)
        end select
      end do
      end do
      end do
      return
    end if

!$omp target update from(BsYstt, BsYend)

    if (n2m /= MPI_PROC_NULL) then
      call MPI_IRECV(BrYstt, size(BrYstt), MPI_DOUBLE_PRECISION, n2m, 2100, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
      call MPI_ISEND(BsYstt, size(BsYstt), MPI_DOUBLE_PRECISION, n2m, 2200, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
    else
      if (boundary_yin == reflection) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,mgn; do i=1,in-1
          BrYstt(i,j,k,1:nbc) = BsYstt(i,mgn-j+1,k,1:nbc)
          BrYstt(i,j,k,nv2) = -BrYstt(i,j,k,nv2)
        end do; end do; end do
      else if (boundary_yin == outflow) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,mgn; do i=1,in-1
          BrYstt(i,j,k,1:nbc) = BsYstt(i,1,k,1:nbc)
        end do; end do; end do
      end if
    end if

    if (n2p /= MPI_PROC_NULL) then
      call MPI_IRECV(BrYend, size(BrYend), MPI_DOUBLE_PRECISION, n2p, 2200, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
      call MPI_ISEND(BsYend, size(BsYend), MPI_DOUBLE_PRECISION, n2p, 2100, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
    else
      if (boundary_yout == reflection) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,mgn; do i=1,in-1
          BrYend(i,j,k,1:nbc) = BsYend(i,mgn-j+1,k,1:nbc)
          BrYend(i,j,k,nv2) = -BrYend(i,j,k,nv2)
        end do; end do; end do
      else if (boundary_yout == outflow) then
!$omp target teams distribute parallel do collapse(3)
        do k=1,kn-1; do j=1,mgn; do i=1,in-1
          BrYend(i,j,k,1:nbc) = BsYend(i,mgn,k,1:nbc)
        end do; end do; end do
      end if
    end if

    if (nreq /= 0) then
      call MPI_WAITALL(nreq, req, stat, ierr)
      if (n2m /= MPI_PROC_NULL) then
!$omp target update to(BrYstt)
      end if
      if (n2p /= MPI_PROC_NULL) then
!$omp target update to(BrYend)
      end if
      nreq = 0
    end if

  end subroutine YbcSendRecv


  subroutine ZbcSendRecv()
    use mpimod
    implicit none
    integer :: i,j,k

    if (ntiles(3) == 1) then
!$omp target teams distribute parallel do collapse(3)
      do j=1,jn-1
      do i=1,in-1
      do k=1,mgn
        select case(boundary_zin)
        case(periodicb)
          BrZstt(i,j,k,1:nbc) = BsZend(i,j,k,1:nbc)
        case(reflection)
          BrZstt(i,j,k,1:nbc) = BsZstt(i,j,mgn-k+1,1:nbc)
          BrZstt(i,j,k,nv3) = -BrZstt(i,j,k,nv3)
        case(outflow)
          BrZstt(i,j,k,1:nbc) = BsZstt(i,j,1,1:nbc)
        end select

        select case(boundary_zout)
        case(periodicb)
          BrZend(i,j,k,1:nbc) = BsZstt(i,j,k,1:nbc)
        case(reflection)
          BrZend(i,j,k,1:nbc) = BsZend(i,j,mgn-k+1,1:nbc)
          BrZend(i,j,k,nv3) = -BrZend(i,j,k,nv3)
        case(outflow)
          BrZend(i,j,k,1:nbc) = BsZend(i,j,mgn,1:nbc)
        end select
      end do
      end do
      end do
      return
    end if

!$omp target update from(BsZstt, BsZend)

    if (n3m /= MPI_PROC_NULL) then
      call MPI_IRECV(BrZstt, size(BrZstt), MPI_DOUBLE_PRECISION, n3m, 3100, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
      call MPI_ISEND(BsZstt, size(BsZstt), MPI_DOUBLE_PRECISION, n3m, 3200, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
    else
      if (boundary_zin == reflection) then
!$omp target teams distribute parallel do collapse(3)
        do j=1,jn-1; do i=1,in-1; do k=1,mgn
          BrZstt(i,j,k,1:nbc) = BsZstt(i,j,mgn-k+1,1:nbc)
          BrZstt(i,j,k,nv3) = -BrZstt(i,j,k,nv3)
        end do; end do; end do
      else if (boundary_zin == outflow) then
!$omp target teams distribute parallel do collapse(3)
        do j=1,jn-1; do i=1,in-1; do k=1,mgn
          BrZstt(i,j,k,1:nbc) = BsZstt(i,j,1,1:nbc)
        end do; end do; end do
      end if
    end if

    if (n3p /= MPI_PROC_NULL) then
      call MPI_IRECV(BrZend, size(BrZend), MPI_DOUBLE_PRECISION, n3p, 3200, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
      call MPI_ISEND(BsZend, size(BsZend), MPI_DOUBLE_PRECISION, n3p, 3100, comm3d, req(nreq+1), ierr)
      nreq = nreq + 1
    else
      if (boundary_zout == reflection) then
!$omp target teams distribute parallel do collapse(3)
        do j=1,jn-1; do i=1,in-1; do k=1,mgn
          BrZend(i,j,k,1:nbc) = BsZend(i,j,mgn-k+1,1:nbc)
          BrZend(i,j,k,nv3) = -BrZend(i,j,k,nv3)
        end do; end do; end do
      else if (boundary_zout == outflow) then
!$omp target teams distribute parallel do collapse(3)
        do j=1,jn-1; do i=1,in-1; do k=1,mgn
          BrZend(i,j,k,1:nbc) = BsZend(i,j,mgn,1:nbc)
        end do; end do; end do
      end if
    end if

    if (nreq /= 0) then
      call MPI_WAITALL(nreq, req, stat, ierr)
      if (n3m /= MPI_PROC_NULL) then
!$omp target update to(BrZstt)
      end if
      if (n3p /= MPI_PROC_NULL) then
!$omp target update to(BrZend)
      end if
      nreq = 0
    end if

  end subroutine ZbcSendRecv

end module boundarymod
