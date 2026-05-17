module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1b,x2b,x3b
    real(8),dimension(:),allocatable:: x1a,x2a,x3a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p
    real(8),dimension(:,:,:),allocatable:: b1,b2,b3,bp
    real(8),dimension(:,:,:),allocatable:: vor1,vor2,vor3
    real(8),dimension(:,:,:),allocatable:: jcd1,jcd2,jcd3
    real(8),dimension(:,:,:),allocatable:: kin,hk
    real(8),dimension(:,:,:),allocatable:: mag,hmm,hcr

    real(8):: dx,dy,dz
end module fieldmod

program data_analysis
  use fieldmod
  implicit none
  integer:: fbeg, fend
  logical:: flag
  integer,parameter:: unitcon=100

  INQUIRE(FILE ="control.dat",EXIST = flag)
  if(flag) then
     open (unitcon,file="control.dat" &
     &        ,status='old',form='formatted')
     read (unitcon,*) fbeg,fend
     close(unitcon)
  endif

  FILENUMBER: do incr  = fbeg,fend
     write(6,*) "file number",incr
     call ReadData
     call Vorticity
     call Snap2D
     call VISITOUT3D
!     call Fourier
!     call Probability
  enddo FILENUMBER

  stop
end program data_analysis


subroutine ReadData
  !! Read binary data.
  !! Grid data contains the information of cell edge.
  !! Variable data contains the information of cell center.
  use fieldmod
  implicit none   
  character(20),parameter::dirname="./bindata/"
  character(40)::filename
  integer::unitunf,unitvar,unitgrd
  integer,parameter:: gs=1
  character(8)::dummy
  logical,save:: is_inited
  real(8),dimension(:,:,:,:),allocatable,save:: data3D
  real(8),dimension(:,:),allocatable,save:: gridX,gridY,gridZ
  data is_inited / .false. /
  integer:: i,j,k

  write(filename,'(a3,i5.5,a4)')"unf",incr,".dat"
  filename = trim(dirname)//filename
  open(newunit=unitunf,file=filename,status='old',form='formatted',action="READ")
  read(unitunf,*) dummy,time,dt
  read(unitunf,*) dummy,izone
  read(unitunf,*) dummy,jzone
  read(unitunf,*) dummy,kzone
  close(unitunf)
  
  print *, "Reading data T=",time,"dt",dt
  print *, "grid",izone,jzone,kzone
 
  in = izone + 2*gs !! +2*gs corresponds to  the inner and outer ghost grid
  jn = jzone + 2*gs !! +2*gs corresponds to  the inner and outer ghost grid
  kn = kzone + 2*gs !! +2*gs corresponds to  the inner and outer ghost grid

  is = 1 + gs 
  js = 1 + gs
  ks = 1 + gs 

  ie = in-gs
  je = jn-gs
  ke = kn-gs

  !> Read grid data.
  if(.not. is_inited)then
     allocate( gridX(2,izone+1))
     allocate( gridY(2,jzone+1))
     allocate( gridZ(2,kzone+1))
     
     filename = trim(dirname)//"g1dDT"
     print *, filename
     open(newunit=unitgrd,file=filename,status='old',form='unformatted',access="stream",action="READ")
     read(unitgrd) gridX(:,:)
     close(unitgrd)
     
     filename = trim(dirname)//"g2dDT"
     open(newunit=unitgrd,file=filename,status='old',form='unformatted',access="stream",action="READ")
     read(unitgrd) gridY(:,:)
     close(unitgrd)
     
     filename = trim(dirname)//"g3dDT"
     open(newunit=unitgrd,file=filename,status='old',form='unformatted',access="stream",action="READ")
     read(unitgrd) gridZ(:,:)
     close(unitgrd)
  
  endif

  !> Map grid data.
  if(.not. is_inited)then
     allocate( x1b(in),x1a(in+1))
     allocate( x2b(jn),x2a(jn+1))
     allocate( x3b(kn),x3a(kn+1))
     
     x1b(is:ie  ) = gridX(1,1:izone  )
     x1a(is:ie+1) = gridX(2,1:izone+1)
     
     x2b(js:je  ) = gridY(1,1:jzone  )
     x2a(js:je+1) = gridY(2,1:jzone+1)
     
     x3b(ks:ke  ) = gridZ(1,1:kzone  )
     x3a(ks:ke+1) = gridZ(2,1:kzone+1)
  endif
  
  
  !> Read variable data.
  if(.not. is_inited)then
     allocate( data3D(9,izone,jzone,kzone))
  endif
  
  write(filename,'(a6,i5.5)')"d3dDT.",incr
  filename = trim(dirname)//filename
  open(newunit=unitvar,file=filename,status='old',form='unformatted', access='stream')
  read(unitvar) data3D(1:9,1:izone,1:jzone,1:kzone)
  close(unitvar)
  
   if(.not. is_inited)then      
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate(b1(in,jn,kn))
     allocate(b2(in,jn,kn))
     allocate(b3(in,jn,kn))
     allocate(bp(in,jn,kn))
     allocate( p(in,jn,kn))
  endif

   d(is:ie,js:je,ks:ke) = data3D(1,1:izone,1:jzone,1:kzone)
  v1(is:ie,js:je,ks:ke) = data3D(2,1:izone,1:jzone,1:kzone)
  v2(is:ie,js:je,ks:ke) = data3D(3,1:izone,1:jzone,1:kzone)
  v3(is:ie,js:je,ks:ke) = data3D(4,1:izone,1:jzone,1:kzone)
  b1(is:ie,js:je,ks:ke) = data3D(5,1:izone,1:jzone,1:kzone)
  b2(is:ie,js:je,ks:ke) = data3D(6,1:izone,1:jzone,1:kzone)
  b3(is:ie,js:je,ks:ke) = data3D(7,1:izone,1:jzone,1:kzone)
  bp(is:ie,js:je,ks:ke) = data3D(8,1:izone,1:jzone,1:kzone)
   p(is:ie,js:je,ks:ke) = data3D(9,1:izone,1:jzone,1:kzone)

  !> Boundary condition
   if(.not. is_inited)then
      dx = x1a(is+1)-x1a(is)
      dy = x2a(js+1)-x2a(js)
      dz = x3a(ks+1)-x3a(ks)
      do i=1,gs
         x1b(is-i) = x1b(is) - i *dx
         x1a(is-i) = x1a(is) - i *dx
      enddo
      do j=1,gs
         x2b(js-j) = x2b(js) - j *dy
         x2a(js-j) = x2a(js) - j *dy
      enddo
      do k=1,gs
         x3b(ks-k) = x3b(ks) - k *dz
         x3a(ks-k) = x3a(ks) - k *dz
      enddo

      do i=1,gs
         x1b(ie  +i) = x1b(ie  ) + i *dx
         x1a(ie+1+i) = x1a(ie+1) + i *dx
      enddo
      do j=1,gs
         x2b(je  +j) = x2b(je  ) + j *dy
         x2a(je+1+j) = x2a(je+1) + j *dy
      enddo
      do k=1,gs
         x3b(ke  +k) = x3b(ke  ) + k *dz
         x3a(ke+1+k) = x3a(ke+1) + k *dz
      enddo
   endif

   !> periodic boundary
   do i=1,gs
       d(is-i,js:je,ks:ke) =  d(ie+1-i,js:je,ks:ke)
      v1(is-i,js:je,ks:ke) = v1(ie+1-i,js:je,ks:ke)
      v2(is-i,js:je,ks:ke) = v2(ie+1-i,js:je,ks:ke)
      v3(is-i,js:je,ks:ke) = v3(ie+1-i,js:je,ks:ke)
      b1(is-i,js:je,ks:ke) = b1(ie+1-i,js:je,ks:ke)
      b2(is-i,js:je,ks:ke) = b2(ie+1-i,js:je,ks:ke)
      b3(is-i,js:je,ks:ke) = b3(ie+1-i,js:je,ks:ke)
      bp(is-i,js:je,ks:ke) = bp(ie+1-i,js:je,ks:ke)
       p(is-i,js:je,ks:ke) =  p(ie+1-i,js:je,ks:ke)
    
       d(ie+i,js:je,ks:ke) =  d(is-1+i,js:je,ks:ke)
      v1(ie+i,js:je,ks:ke) = v1(is-1+i,js:je,ks:ke)
      v2(ie+i,js:je,ks:ke) = v2(is-1+i,js:je,ks:ke)
      v3(ie+i,js:je,ks:ke) = v3(is-1+i,js:je,ks:ke)
      b1(ie+i,js:je,ks:ke) = b1(is-1+i,js:je,ks:ke)
      b2(ie+i,js:je,ks:ke) = b2(is-1+i,js:je,ks:ke)
      b3(ie+i,js:je,ks:ke) = b3(is-1+i,js:je,ks:ke)
      bp(ie+i,js:je,ks:ke) = bp(is-1+i,js:je,ks:ke)
       p(ie+i,js:je,ks:ke) =  p(is-1+i,js:je,ks:ke)
   enddo


     do j=1,gs
       d(is:ie,js-j,ks:ke) =  d(is:ie,je+1-j,ks:ke)
      v1(is:ie,js-j,ks:ke) = v1(is:ie,je+1-j,ks:ke)
      v2(is:ie,js-j,ks:ke) = v2(is:ie,je+1-j,ks:ke)
      v3(is:ie,js-j,ks:ke) = v3(is:ie,je+1-j,ks:ke)
      b1(is:ie,js-j,ks:ke) = b1(is:ie,je+1-j,ks:ke)
      b2(is:ie,js-j,ks:ke) = b2(is:ie,je+1-j,ks:ke)
      b3(is:ie,js-j,ks:ke) = b3(is:ie,je+1-j,ks:ke)
      bp(is:ie,js-j,ks:ke) = bp(is:ie,je+1-j,ks:ke)
       p(is:ie,js-j,ks:ke) =  p(is:ie,je+1-j,ks:ke)
    
       d(is:ie,je+j,ks:ke) =  d(is:ie,js-1+j,ks:ke)
      v1(is:ie,je+j,ks:ke) = v1(is:ie,js-1+j,ks:ke)
      v2(is:ie,je+j,ks:ke) = v2(is:ie,js-1+j,ks:ke)
      v3(is:ie,je+j,ks:ke) = v3(is:ie,js-1+j,ks:ke)
      b1(is:ie,je+j,ks:ke) = b1(is:ie,js-1+j,ks:ke)
      b2(is:ie,je+j,ks:ke) = b2(is:ie,js-1+j,ks:ke)
      b3(is:ie,je+j,ks:ke) = b3(is:ie,js-1+j,ks:ke)
      bp(is:ie,je+j,ks:ke) = bp(is:ie,js-1+j,ks:ke)
       p(is:ie,je+j,ks:ke) =  p(is:ie,js-1+j,ks:ke)
   enddo

   do k=1,gs
       d(is:ie,js:je,ks-k) =  d(is:ie,js:je,ke+1-k)
      v1(is:ie,js:je,ks-k) = v1(is:ie,js:je,ke+1-k)
      v2(is:ie,js:je,ks-k) = v2(is:ie,js:je,ke+1-k)
      v3(is:ie,js:je,ks-k) = v3(is:ie,js:je,ke+1-k)
      b1(is:ie,js:je,ks-k) = b1(is:ie,js:je,ke+1-k)
      b2(is:ie,js:je,ks-k) = b2(is:ie,js:je,ke+1-k)
      b3(is:ie,js:je,ks-k) = b3(is:ie,js:je,ke+1-k)
      bp(is:ie,js:je,ks-k) = bp(is:ie,js:je,ke+1-k)
       p(is:ie,js:je,ks-k) =  p(is:ie,js:je,ke+1-k)
    
       d(is:ie,js:je,ke+k) =  d(is:ie,js:je,ks-1+k)
      v1(is:ie,js:je,ke+k) = v1(is:ie,js:je,ks-1+k)
      v2(is:ie,js:je,ke+k) = v2(is:ie,js:je,ks-1+k)
      v3(is:ie,js:je,ke+k) = v3(is:ie,js:je,ks-1+k)
      b1(is:ie,js:je,ke+k) = b1(is:ie,js:je,ks-1+k)
      b2(is:ie,js:je,ke+k) = b2(is:ie,js:je,ks-1+k)
      b3(is:ie,js:je,ke+k) = b3(is:ie,js:je,ks-1+k)
      bp(is:ie,js:je,ke+k) = bp(is:ie,js:je,ks-1+k)
       p(is:ie,js:je,ke+k) =  p(is:ie,js:je,ks-1+k)
    enddo
 
  is_inited = .true.
    
  return
end subroutine ReadData

subroutine Vorticity
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     allocate( vor1(in,jn,kn))
     allocate( vor2(in,jn,kn))
     allocate( vor3(in,jn,kn))
     allocate( jcd1(in,jn,kn))
     allocate( jcd2(in,jn,kn))
     allocate( jcd3(in,jn,kn))
     allocate(  kin(in,jn,kn))
     allocate(   hk(in,jn,kn))
     allocate(  mag(in,jn,kn))
     allocate(  hmm(in,jn,kn))
     allocate(  hcr(in,jn,kn))
     is_inited = .true.
  endif

  do k=ks,ke
  do j=js,je
  do i=is,ie
     vor1(i,j,k)= (v3(i,j  ,k)-v3(i,j-1,k))/dy*0.5 &
                &+(v3(i,j+1,k)-v3(i,j  ,k))/dy*0.5 &
                &-(v2(i,j,k  )-v2(i,j,k-1))/dz*0.5 &
                &-(v2(i,j,k+1)-v2(i,j,k  ))/dz*0.5 
     vor2(i,j,k)= (v1(i,j,k  )-v1(i,j,k-1))/dz*0.5 &
                &+(v1(i,j,k+1)-v1(i,j,k  ))/dz*0.5 &
                &-(v3(i  ,j,k)-v3(i-1,j,k))/dx*0.5 &
                &-(v3(i+1,j,k)-v3(i  ,j,k))/dx*0.5 
     vor3(i,j,k)= (v2(i  ,j,k)-v2(i-1,j,k))/dx*0.5 &
                &+(v2(i+1,j,k)-v2(i  ,j,k))/dx*0.5 &
                &-(v1(i,j  ,k)-v1(i,j-1,k))/dy*0.5 &
                &-(v1(i,j+1,k)-v1(i,j  ,k))/dy*0.5

      kin(i,j,k)= 0.5d0*d(i,j,k)*( v1(i,j,k)*v1(i,j,k) &
                &                 +v2(i,j,k)*v2(i,j,k) &
                &                 +v3(i,j,k)*v3(i,j,k))

       hk(i,j,k)=  v1(i,j,k)*vor1(i,j,k) &
                & +v2(i,j,k)*vor2(i,j,k) &
                & +v3(i,j,k)*vor3(i,j,k)

     jcd1(i,j,k)= (b3(i,j  ,k)-b3(i,j-1,k))/dy*0.5 &
                &+(b3(i,j+1,k)-b3(i,j  ,k))/dy*0.5 &
                &-(b2(i,j,k  )-b2(i,j,k-1))/dz*0.5 &
                &-(b2(i,j,k+1)-b2(i,j,k  ))/dz*0.5 
     jcd2(i,j,k)= (b1(i,j,k  )-b1(i,j,k-1))/dz*0.5 &
                &+(b1(i,j,k+1)-b1(i,j,k  ))/dz*0.5 &
                &-(b3(i  ,j,k)-b3(i-1,j,k))/dx*0.5 &
                &-(b3(i+1,j,k)-b3(i  ,j,k))/dx*0.5 
     jcd3(i,j,k)= (b2(i  ,j,k)-b2(i-1,j,k))/dx*0.5 &
                &+(b2(i+1,j,k)-b2(i  ,j,k))/dx*0.5 &
                &-(b1(i,j  ,k)-b1(i,j-1,k))/dy*0.5 &
                &-(b1(i,j+1,k)-b1(i,j  ,k))/dy*0.5

      mag(i,j,k)= 0.5d0*( b1(i,j,k)*b1(i,j,k) &
                &        +b2(i,j,k)*b2(i,j,k) &
                &        +b3(i,j,k)*b3(i,j,k))

       hmm(i,j,k)=  b1(i,j,k)*jcd1(i,j,k) &
                 & +b2(i,j,k)*jcd2(i,j,k) &
                 & +b3(i,j,k)*jcd3(i,j,k)
      hcr(i,j,k)=    ( b1(i,j,k)*v1(i,j,k) &
                &     +b2(i,j,k)*v2(i,j,k) &
                &     +b3(i,j,k)*v3(i,j,k))



  enddo
  enddo
  enddo

  return
end subroutine Vorticity

      subroutine VISITOUT3D
      use fieldmod
      use hdf5
      implicit none 
      integer,parameter::ndim=3                   ! max size of the dimension for output file
      integer(hid_t) :: file_id                   ! file identifier
      integer ::   rank                           ! dataset rank
      integer(hsize_t), dimension(ndim) :: dims   ! dataset dimensions
      character(80) :: dsetname                   ! dataset name
      integer :: error
      character*(80) :: fname,fpath
      character*(80) :: fnxmf,fpxmf
      integer unitxmf
      integer,save:: imax,jmax,kmax
      character*(80) :: DIM3D,DIM1D,TIMEEXP

      character(40),parameter:: vstdir = "./vstdata/" 

      character(2):: id

      integer,parameter:: nfld=5
      integer,save:: vnum
      data vnum /0/

! Variable to esatimate NODE value

      real*8, dimension(:,:,:,:), allocatable,save:: V3DCEND
      real*4, dimension(:,:,:,:), allocatable,save:: V3DCENF
      real*4,dimension(:),allocatable,save::XF,YF,ZF

      integer:: i,j,k

      logical,save::is_inited
      data is_inited / .false. /

      id="tb"

! Initialize
      if(.not. is_inited)then
         imax = in-2*igs
         jmax = jn-2*jgs
         kmax = kn-2*kgs

         allocate(XF(imax+1))
         allocate(YF(jmax+1))
         allocate(ZF(kmax+1))
         allocate(V3DCENF   (1:imax  ,1:jmax  ,1:kmax  ,nfld))
         
         XF(1:imax+1)= real(x1a(1:imax+1))
         YF(1:jmax+1)= real(x2a(1:jmax+1))
         ZF(1:kmax+1)= real(x3a(1:kmax+1))
  
         call makedirs(vstdir)

         is_inited = .true.
      endif
      
!      if(.not. mod(incr,10) ==0 )return
      vnum=vnum+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing HDF5 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine the output filename
            write(fname,"(a2,i5.5,a3)") &
     &  id,vnum,'.h5'
      fpath =  trim(vstdir)//fname
! Open HDF5
      CALL h5open_f (error)
! Note:
! H5F_ACC_TRUNC_F: delete the file if that is exist
! H5F_ACC_EXCL_F : fail writing to the file if that  is exist
      CALL h5fcreate_f(fpath, H5F_ACC_TRUNC_F, file_id, error)
      dims(:) = 0
      rank = 1
      dims(1) = imax+1
      dsetname="/x"
      call  hdf5OUT(file_id,rank,dims,dsetname,XF)
      rank=1
      dims(1) = jmax+1
      dsetname="/y"
      call  hdf5OUT(file_id,rank,dims,dsetname,YF)
      rank=1
      dims(1) = kmax+1
      dsetname="/z"
      call  hdf5OUT(file_id,rank,dims,dsetname,ZF)

      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         V3DCENF(i,j,k,1) = real( kin(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,2) = real(  hk(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,3) = real( mag(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,4) = real( hmm(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,5) = real( hcr(igs+i,jgs+j,kgs+k) )
      enddo
      enddo
      enddo

      rank=3
      dims(1) = imax
      dims(2) = jmax
      dims(3) = kmax
      dsetname="/kinene"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,1))
      dsetname="/kinhel"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,2))
      dsetname="/magene"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,3))
      dsetname="/maghlm"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,4))
      dsetname="/crshel"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,5))

! Close the File
      call  h5Fclose_f(file_id,error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing XDMF file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine the output filename
            write(fnxmf,"(a2,i5.5,a4)") & 
     &  id,vnum,'.xmf'
      fpxmf =  trim(vstdir)//fnxmf
! FileOpen
      unitxmf=1221
      open(unitxmf,file=fpxmf,status='unknown',form='formatted')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Header
      write(unitxmf,'(a)')'<?xml version="1.0" ?>'
      write(unitxmf,'(a)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(unitxmf,'(a)')'<Xdmf Version="2.0">'
      write(unitxmf,'(a)')'  <Domain>'
      write(unitxmf,'(a)')'    <Grid Name="mesh" GridType="Uniform">'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time
      write(TIMEEXP,"(I0)") int((time))
      write(unitxmf,'(a)')'      <Time Value="'//trim(TIMEEXP)//'"/>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coordinate
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax+1,jmax+1,imax+1
      write(unitxmf,'(a)')'      <Topology' &
     &               //' TopologyType="3DRectMesh" NumberOfElements="' &
     &               //trim(DIM3D)//'"/>'
      write(unitxmf,'(a)')'      <Geometry GeometryType="VXVYVZ">'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x coordinate
      write(DIM1D,"(I0)") imax+1
      dsetname = "x"
      write(unitxmf,'(a)')&
     &                '        <DataItem Dimensions="' // trim(DIM1D)&
     &              //'" Name="' // trim(dsetname) // '"'
      write(unitxmf,'(a)') '          ' // &
     &               'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  ' // trim(fname) //":/x"
      write(unitxmf,'(a)')'        </DataItem>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y coordinate
      write(DIM1D,"(I0)") jmax+1
      dsetname="y"
      write(unitxmf,*)'        <DataItem Dimensions="'//trim(DIM1D)&
     &              //'" Name="' // trim(dsetname) // '"'
      write(unitxmf,*) '          ' // &
     &              'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/y"
      write(unitxmf,'(a)')'        </DataItem>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z coordinate
      write(DIM1D,"(I0)") kmax+1
      dsetname="z"
      write(unitxmf,'(a)')'        <DataItem Dimensions="'//trim(DIM1D) &
     &              //'" Name="'//trim(dsetname) // '"'
      write(unitxmf,'(a)') '          ' // &
     &              'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/z"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Geometry>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kinetic energy
      dsetname="kinene"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/kinene"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kinematic helicity
      dsetname="kinhel"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/kinhel"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! magnetic energy
      dsetname="magene"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/magene"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! magnetic helicity mimic
      dsetname="maghlm"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/maghlm"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cross helicity
      dsetname="crshel"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/crshel"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Footer
      write(unitxmf,'(a)')'    </Grid>'
      write(unitxmf,'(a)')'  </Domain>'

      write(unitxmf,'(a)')'</Xdmf>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      close(unitxmf)
      return
      end subroutine VISITOUT3D

      subroutine hdf5OUT(file_id,rank,dims,dsetname,dset)
      use hdf5
      implicit none
      integer(hid_t),intent(in) :: file_id                    ! file identifier
      integer,intent(in)  :: rank                             ! dataset rank
      integer(hsize_t), dimension(rank),intent(in)  :: dims   ! dataset dimensions
      character(80),intent(in):: dsetname                     ! dataset name
      real(4),intent(in) :: dset
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
! Create the data space for the data set
      call h5Screate_simple_f(rank, dims, dspace_id, error)
! Get dset_id for data set
      call h5Dcreate_f(file_id,dsetname,h5t_native_real,dspace_id, &
     &                 dset_id,error)
! Write the data
      call h5Dwrite_f(dset_id, h5t_native_real, dset, dims, error)
! Close the data id
      call h5Dclose_f(dset_id, error)
! Close the data space
      call h5Sclose_f(dspace_id, error)

      return
      end subroutine hdf5out

      subroutine makedirs(outdir)
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs

subroutine Snap2D
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /

  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')
  k=ks
  write(unitvor,*) "# ",time
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(6(1x,E12.3))') x1b(i),x2b(j),vor3(i,j,k),jcd3(i,j,k),kin(i,j,k),mag(i,j,k)
  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Snap2D
 
