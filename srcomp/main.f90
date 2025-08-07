program main
  use omp_lib
  use basicmod
  use mpimod
  use boundarymod
  implicit none
  real(8)::time_begin,time_end
  logical::is_final
  logical,parameter::nooutput=.true.
  logical,parameter::debugmode=.true.
  data is_final /.false./
  call InitializeMPI
  if(myid_w == 0) print *, "setup grids and fields"
  if(myid_w == 0) print *, "grid size for x y z",ngrid1*ntiles(1),ngrid2*ntiles(2),ngrid3*ntiles(3)
  if(myid_w == 0 .and. nooutput ) print *, "Intermediate results are not outputed"
  call GenerateGrid
  call GenerateProblem
  call ConsvVariable
  if(myid_w == 0) print *, "entering main loop"
! main loop
  if(myid_w == 0 .and. .not. nooutput )                        print *,"step ","time ","dt"
  time_begin = omp_get_wtime()
  mloop: do nhy=1,nhymax
     
     if(debugmode) print *,"TimestepControl"
     call TimestepControl
     if(mod(nhy,nhydis) .eq. 0  .and. .not. nooutput .and. myid_w == 0) print *,nhy,time,dt
     if(debugmode) print *,"BoundaryCondition"
     call BoundaryCondition
     call StateVevtor
     if(debugmode) print *,"EvaulateCh"
     call EvaulateCh
     if(debugmode) print *,"NumericalFlux"
     call NumericalFlux1
     call NumericalFlux2
     call NumericalFlux3
     if(debugmode) print *,"UpdateConsv"
     call UpdateConsv
     if(debugmode) print *,"DampPsi"
     call DampPsi
     if(debugmode) print *,"PrimVariable"
     call PrimVariable
     time=time+dt
     if(.not. nooutput ) call Output(.false.)
     if(time > timemax) exit mloop
  enddo mloop

  time_end = omp_get_wtime()
      
  if(myid_w == 0) print *, "sim time [s]:", time_end-time_begin
  if(myid_w == 0) print *, "time/count/cell", (time_end-time_begin)/(ngrid1*ngrid2*ngrid3)/nhymax
  
  is_final = .true.
  call Output(.true.)

  call FinalizeMPI
  if(myid_w == 0) print *, "program has been finished"
  
end program main


subroutine GenerateGrid
  use basicmod
  use mpimod
  implicit none
  real(8)::dx,dy,dz
  real(8)::x1minloc,x1maxloc
  real(8)::x2minloc,x2maxloc
  real(8)::x3minloc,x3maxloc
  integer::i,j,k
  
  ! x coordinates
      
  x1minloc = x1min + (x1max-x1min)/ntiles(1)* coords(1)
  x1maxloc = x1min + (x1max-x1min)/ntiles(1)*(coords(1)+1)
  
  dx=(x1maxloc-x1minloc)/dble(ngrid1)
  do i=1,in
     x1a(i) = dx*(i-(mgn+1))+x1minloc
  enddo
  do i=1,in-1
     x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
  enddo
 
  ! y coordinates
  x2minloc = x2min + (x2max-x2min)/ntiles(2)* coords(2)
  x2maxloc = x2min + (x2max-x2min)/ntiles(2)*(coords(2)+1)
 
  dy=(x2maxloc-x2minloc)/dble(ngrid2)
  do j=1,jn
     x2a(j) = dy*(j-(mgn+1))+x2minloc
  enddo

  do j=1,jn-1
     x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
  enddo
  
  ! z coordinates
  x3minloc = x3min + (x3max-x3min)/ntiles(3)* coords(3)
  x3maxloc = x3min + (x3max-x3min)/ntiles(3)*(coords(3)+1)
 
  dz=(x3maxloc-x3minloc)/ngrid3
  do k=1,kn
     x3a(k) = dz*(k-(mgn+1))+x3minloc
  enddo
  do k=1,kn-1
     x3b(k) = 0.5d0*(x3a(k+1)+x3a(k))
  enddo

!$omp target update to(x1a,x1b)
!$omp target update to(x2a,x2b)
!$omp target update to(x3a,x3b)
  
  return
end subroutine GenerateGrid

subroutine GenerateProblem
  use basicmod
  use eosmod
  use mpimod
  use boundarymod
  implicit none
  integer::i,j,k

  real(8),parameter::pi=acos(-1.0d0)

  real(8)::Ahl,Bhl,Chl
  real(8),parameter::k_ini=2.0d0
  
  real(8),parameter:: ekin = 2.0d0
  real(8),parameter:: emag = 2.0d0
  real(8),parameter:: eint = 1.0d0
  real(8),parameter:: d0 = 1.0d0
  real(8),parameter:: v0 = sqrt(ekin*2.d0/d0)
  real(8),parameter:: b0 = sqrt(emag*2.0)
  real(8)          :: p0
  real(8),parameter:: eps = 1.0d-1
  real(8),parameter:: deltax = 0.1d0,deltay = 0.2d0,deltaz = 0.3d0 ! randam phase

  integer::seedsize
  integer,allocatable:: seed(:)
  real(8)::x

  call random_seed(size=seedsize)
!  print *,"seed size",seedsize
  allocate(seed(seedsize))  
  call random_seed(get=seed)

      Ahl = 0.5d0
      Bhl = 0.5d0
      Chl = 0.5d0

      d(:,:,:) = d0
! adiabatic
!       p0= eint/(gam-1.0d0)
! isotermal
       csiso= sqrt(eint/d0)
       p0 = d0 *csiso**2       
!$omp target update to(csiso)
       
      do k=ks,ke
      do j=js,je
      do i=is,ie
         v1(i,j,k) = v0*(  Ahl*sin(2.0d0*pi*(k_ini*x3b(k)/(x3max-x3min)+deltaz)) &
   &                     + Chl*cos(2.0d0*pi*(k_ini*x2b(j)/(x2max-x2min)+deltay)))
         v2(i,j,k) = v0*(  Bhl*sin(2.0d0*pi*(k_ini*x1b(i)/(x1max-x1min)+deltax)) &
   &                     + Ahl*cos(2.0d0*pi*(k_ini*x3b(k)/(x3max-x3min)+deltaz)))
         v3(i,j,k) = v0*(  Chl*sin(2.0d0*pi*(k_ini*x2b(j)/(x2max-x2min)+deltay)) &
   &                     + Bhl*cos(2.0d0*pi*(k_ini*x1b(i)/(x1max-x1min)+deltax)))

          p(i,j,k) = p0

         b1(i,j,k) = b0*(  Ahl*sin(2.0d0*pi*(k_ini*x3b(k)/(x3max-x3min))) &
   &                     + Chl*cos(2.0d0*pi*(k_ini*x2b(j)/(x2max-x2min))))
         b2(i,j,k) = b0*(  Bhl*sin(2.0d0*pi*(k_ini*x1b(i)/(x1max-x1min))) &
   &                     + Ahl*cos(2.0d0*pi*(k_ini*x3b(k)/(x3max-x3min))))
         b3(i,j,k) = b0*(  Chl*sin(2.0d0*pi*(k_ini*x2b(j)/(x2max-x2min))) &
   &                     + Bhl*cos(2.0d0*pi*(k_ini*x1b(i)/(x1max-x1min))))

         call random_number(x)
         v1(i,j,k) = v1(i,j,k)*(1.0d0+eps*(x-0.5d0))
         call random_number(x)
         v2(i,j,k) = v2(i,j,k)*(1.0d0+eps*(x-0.5d0))
         call random_number(x)
         v3(i,j,k) = v3(i,j,k)*(1.0d0+eps*(x-0.5d0))
      enddo
      enddo
      enddo


      do k=ks,ke
      do j=js,je
      do i=is,ie
! adiabatic
!          ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
!          cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
! isotermal
          ei(i,j,k) = p(i,j,k)
          cs(i,j,k) = csiso
      enddo
      enddo
      enddo
      
      if(myid_w ==0 )print *,"initial profile is set"
      call BoundaryCondition

!$omp target update to(d,v1,v2,v3)
!$omp target update to(p,ei,cs)
!$omp target update to(b1,b2,b3,bp)
      
  return
end subroutine GenerateProblem

! Code was translated using: /gwork0/takiwkkz/MHDturbArxiv/acctoomp/src/intel-application-migration-tool-for-openacc-to-openmp -suppress-openacc main.f90
