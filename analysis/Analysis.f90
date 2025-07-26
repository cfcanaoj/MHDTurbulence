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
    
!$acc declare create(incr)
!$acc declare create(time,dt)
!$acc declare create(in,jn,kn)
!$acc declare create(izone,jzone,kzone)
!$acc declare create(igs,jgs,kgs)
!$acc declare create(is,js,ks)
!$acc declare create(ie,je,ke)
      
!$acc declare create(x1a,x1b)
!$acc declare create(x2a,x2b)
!$acc declare create(x3a,x3b)
      
!$acc declare create(d,v1,v2,v3,p)
!$acc declare create(b1,b2,b3,bp)
!$acc declare create(vor1,vor2,vor3)
!$acc declare create(jcd1,jcd2,jcd3)
!$acc declare create(kin,hk)
!$acc declare create(mag,hmm,hcr)

!$acc declare create(dx,dy,dz)

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

  FILENUMBER: do incr  = fbeg,fend,10
     write(6,*) "file number",incr
     call ReadData
     call Vorticity
     call Fourier
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
!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (x3a,x3b)
!$acc update device (d,v1,v2,v3,p)
!$acc update device (b1,b2,b3,bp)
!$acc update device (dx,dy,dz)

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
!$acc update device (vor1,vor2,vor3)
!$acc update device (jcd1,jcd2,jcd3)
!$acc update device (kin,hk)
!$acc update device (mag,hmm,hcr)
     is_inited = .true.
  endif

!$acc kernels
!$acc loop collapse(3) independent
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
!$acc end kernels

  
!$acc update host (kin,hk)
!$acc update host (mag,hmm,hcr)
  write(6,*)"debug1",kin(is,js,ks),hk(is,js,ks),mag(is,js,ks),hmm(is,js,ks),hcr(is,js,ks)
  return
end subroutine Vorticity

module spctrmod
implicit none
  real(8),dimension(:,:,:,:),allocatable:: X3D
!$acc declare create(X3D)
  integer,parameter:: nk=128
  integer,parameter:: nvar=5
  real(8),dimension(nvar):: Xtot
  real(8),dimension(nk,nk,nk,nvar):: Xhat3DC,Xhat3DS
  real(8),dimension(nk):: kx,ky,kz
  real(8),dimension(nk,nvar):: Xhat1D
  real(8):: kr
  real(8):: dkx,dky,dkz,dkr
  
  real(8) :: pi
!$acc declare create(Xtot)
!$acc declare create(dkx,dky,dkz)
!$acc declare create(kx,ky,kz)
!$acc declare create(Xhat3DC,Xhat3DS,Xhat1D)
!$acc declare create(pi)
end module spctrmod

subroutine Fourier
  use fieldmod
  use spctrmod
  implicit none
  integer::i,j,k,n
  integer::ik,jk,kk,rk
  real(8):: Xtotloc
  real(8):: Xhat3DCloc,Xhat3DSloc
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  integer,parameter::unittot=22
  integer,save::nout
  data nout / 1 /
  logical,save:: is_inited
  data is_inited / .false. /
  
  if(.not. is_inited)then
     allocate( X3D(is:ie,js:je,ks:ke,nvar))
     pi=acos(-1.0d0)
!$acc update device (X3D)
!$acc update device (pi)
     is_inited = .true.
  endif

!$acc kernels
!$acc loop collapse(3) independent
  do k=ks,ke
  do j=js,je
  do i=is,ie
     X3D(i,j,k,1) = kin(i,j,k)
     X3D(i,j,k,2) =  hk(i,j,k)
     X3D(i,j,k,3) = hmm(i,j,k)
     X3D(i,j,k,4) = Hcr(i,j,k)
     X3D(i,j,k,5) = mag(i,j,k)
  enddo
  enddo
  enddo
!$acc end kernels  

!$acc update host (X3D)
  write(6,*)"debug2",X3D(is,js,ks,1:nvar)
 
!$acc kernels
!$acc loop independent private(Xtotloc)
  do n=1,nvar
     Xtotloc = 0.0d0
!$acc loop collapse(3)reduction(+:Xtotloc)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Xtotloc = Xtotloc + X3D(i,j,k,n)*dx*dy*dz
  enddo
  enddo
  enddo
     Xtot(n) = Xtotloc
  enddo
!$acc end kernels
!$acc update host (Xtot)
  
  write(6,*)"debug3",Xtot(1:nvar)

  
  dkx = 1.0d0/(dx*(in-2*igs))
  dky = 1.0d0/(dy*(jn-2*jgs))
  dkz = 1.0d0/(dz*(kn-2*kgs))
!$acc update device (dkx,dky,dkz)
  
  do ik=1,nk
     kx(ik) = ik *dkx
  enddo
  do jk=1,nk
     ky(jk) = jk *dky
  enddo
  do kk=1,nk
     ky(kk) = kk *dkz
  enddo
!$acc update device (kx,ky,kz)

!$acc kernels
!$acc loop collapse(4) independent private(Xhat3DCloc)
  do n=1,nvar
  do kk=1,nk
  do jk=1,nk
  do ik=1,nk 
     Xhat3DCloc = 0.0d0    
!$acc loop collapse(3) reduction(+:Xhat3DCloc)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Xhat3DCloc = Xhat3DCloc &
 &    + X3D(i,j,k,n) &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j)+kz(kk)*x3b(k) )) & 
 &    *dx*dy*dz

  enddo
  enddo
  enddo
     Xhat3DC(ik,jk,kk,n) = Xhat3DCloc
  enddo
  enddo
  enddo
  enddo
!$acc end kernels
!$acc update host (Xhat3DC)

!$acc kernels
!$acc loop collapse(4) independent private(Xhat3DSloc)
  do n=1,nvar  
  do kk=1,nk
  do jk=1,nk
  do ik=1,nk     
     Xhat3DSloc = 0.0d0
!$acc loop collapse(3) reduction(+:Xhat3DSloc)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Xhat3DSloc = Xhat3DSloc &
 &    + X3D(i,j,k,n) &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j)+kz(kk)*x3b(k) )) & 
 &    *dx*dy*dz

  enddo
  enddo
  enddo
     Xhat3DS(ik,jk,kk,n) = Xhat3DSloc 
  enddo
  enddo
  enddo
  enddo
!$acc end kernels
!$acc update host (Xhat3DS)

  Xhat1D(:,:) = 0.0d0
  dkr = dkx/sqrt(3.0d0) ! minimum k
  do kk=1,nk
  do jk=1,nk
  do ik=1,nk
     kr = sqrt(kx(ik)**2 +ky(jk)**2 +kz(kk)**2)
     rk = min(nk,int(kr/dkr))
     Xhat1D(rk,1:nvar) = Xhat1D(rk,1:nvar) + sqrt(   Xhat3DS(ik,jk,kk,1:nvar)**2 & 
                                           &       + Xhat3DC(ik,jk,kk,1:nvar)**2 &
                                           &      )*dkx*dky*dkz
  enddo
  enddo
  enddo

  
  write(filename,'(a3,i5.5,a4)')"spc",nout,".dat"
  filename = trim(dirname)//filename
  open(unitspc,file=filename,status='replace',form='formatted')
  write(unitspc,*) "# ",time
  do rk=1,nk
     write(unitspc,'(7(1x,E12.3))') rk*dkr,Xhat1D(rk,1)/Xtot(1) &
                                         &,Xhat1D(rk,2) &
                                         &,Xhat1D(rk,3) &
                                         &,Xhat1D(rk,4) &
                                         &,Xhat1D(rk,5)/Xtot(5)
  enddo
  close(unitspc)

  write(filename,'(a3,i5.5,a4)')"tot",nout,".dat"
  filename = trim(dirname)//filename
  open(unittot,file=filename,status='replace',form='formatted')
  write(unittot,'(6(1x,E12.3))') time,Xtot(1),Xtot(2),Xtot(3),Xtot(4),Xtot(5)
  close(unittot)

  nout=nout+1
  
  return
end subroutine Fourier
