program main
  use omp_lib
  use basicmod
  use mpimod
  use boundarymod
  implicit none
  real(8)::time_begin,time_end
  logical::is_final
  logical,parameter::nooutput=.true.
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
  if(myid_w == 0 .and. .not. nooutput )                        print *,"step","time","dt"
  time_begin = omp_get_wtime()
  mloop: do nhy=1,nhymax
     call TimestepControl
     if(mod(nhy,nhydis) .eq. 0  .and. .not. nooutput ) print *,nhy,time,dt
     call BoundaryCondition
     call StateVevtor
     call EvaulateCh
     call NumericalFlux1
     call NumericalFlux2
     call NumericalFlux3
     call UpdateConsv
     call DampPsi
     call PrimVariable
     time=time+dt
     if(myid_w == 0 .and. .not. nooutput ) call Output(.false.)
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

!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (x3a,x3b)
  
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
!$acc update device (csiso)
       
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

!$acc update device (d,v1,v2,v3)
!$acc update device (p,ei,cs)
!$acc update device (b1,b2,b3,bp)
      
  return
end subroutine GenerateProblem
      
subroutine Output(is_final)
  !! Output the grid and variables.
  !! The grid data contain the information of cell center and cell edge.
  !! The variable data contains that of  cell center.
  use basicmod
  use mpiiomod
  use mpimod
  implicit none
  integer::i,j,k
  integer::iee,jee,kee
  character(20),parameter::dirname="bindata/"
  character(40)::filename
  real(8),save::tout
  data tout / 0.0d0 /
  integer::nout
  data nout / 1 /
  integer,parameter::unitout=17
  integer,parameter::unitbin=13
  integer,parameter:: gs=1
  integer,parameter:: nvar=9

  logical, intent(in):: is_final

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
     
     nvarg = 2
     nvars = 9
     
     allocate(gridX(nvarg,1:iee-is+1))
     allocate(gridY(nvarg,1:jee-js+1))
     allocate(gridZ(nvarg,1:kee-ks+1))

     allocate(data3D(nvars,ngrid1,ngrid2,ngrid3))
     
     call makedirs("bindata")
     is_inited =.true.
  endif

  if(time .lt. tout+dtout .and. .not. is_final) return

  !> Information of meta data.
  if(myid_w == 0)then
  write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
  filename = trim(dirname)//filename

  open(unitout,file=filename,status='replace',form='formatted')
  write(unitout,'(a2,2(1x,E12.3))') "# ",time,dt
  write(unitout,'(a2,1x,i5)') "# ",ngrid1*ntiles(1)
  write(unitout,'(a2,1x,i5)') "# ",ngrid2*ntiles(2)
  write(unitout,'(a2,1x,i5)') "# ",ngrid3*ntiles(3)
  close(unitout)
  endif

  gridX(1,1:iee-is+1) = x1b(is:iee) !! the final grid point is not necessary but outputed.  
  gridX(2,1:iee-is+1) = x1a(is:iee) !! the final grid is necessary. 
  
  gridY(1,1:jee-js+1) = x2b(js:jee) !! the final grid point is not necessary but outputed.
  gridY(2,1:jee-js+1) = x2a(js:jee) !! the final grid is necessary.

  gridZ(1,1:kee-ks+1) = x3b(ks:kee) !! the final grid point is not necessary but outputed.
  gridZ(2,1:kee-ks+1) = x3a(ks:kee) !! the final grid is necessary.

  !> The cell center value 
  data3D(1,1:npart(1),1:npart(2),1:npart(3)) =  d(is:ie,js:je,ks:ke)
  data3D(2,1:npart(1),1:npart(2),1:npart(3)) = v1(is:ie,js:je,ks:ke)
  data3D(3,1:npart(1),1:npart(2),1:npart(3)) = v2(is:ie,js:je,ks:ke)
  data3D(4,1:npart(1),1:npart(2),1:npart(3)) = v3(is:ie,js:je,ks:ke)
  data3D(5,1:npart(1),1:npart(2),1:npart(3)) = b1(is:ie,js:je,ks:ke)
  data3D(6,1:npart(1),1:npart(2),1:npart(3)) = b2(is:ie,js:je,ks:ke)
  data3D(7,1:npart(1),1:npart(2),1:npart(3)) = b3(is:ie,js:je,ks:ke)
  data3D(8,1:npart(1),1:npart(2),1:npart(3)) = bp(is:ie,js:je,ks:ke)
  data3D(9,1:npart(1),1:npart(2),1:npart(3)) =  p(is:ie,js:je,ks:ke)

  if(myid_w==0)print *, "output:",nout,time

  call MPIOutputBindary(nout)
      
  nout=nout+1
  tout=time
  return
!         write(6,*) "bpf2",nflux2(mbps,i,j,k)
end subroutine Output
      
subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!  write(*, *) trim(command)
  call system(command)
end subroutine makedirs
