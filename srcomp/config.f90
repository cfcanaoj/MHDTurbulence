      module config
      implicit none
      !> Control timescale
      real(8),parameter:: timemax = 15.0d0
      real(8),parameter:: dtout = timemax/100

      !> Control step
      integer,parameter:: nhymax = 600000
      integer,parameter:: nhydis = nhymax/100
      
      !> Control resolution
      integer,parameter:: ngridtotal1 = 150
      integer,parameter:: ngridtotal2 = 150
      integer,parameter:: ngridtotal3 = 150

      !> Control computational region
      real(8),parameter:: x1min = -0.5d0, x1max = 0.5d0
      real(8),parameter:: x2min = -1.0d0, x2max = 1.0d0
      real(8),parameter:: x3min = -0.5d0, x3max = 0.5d0      
      integer,parameter:: ncomp=1 ! kinds of composition
      
      !> Control MPI decomposion
      integer,parameter:: ntiles(3) = [ 2,2,1 ]
      logical,parameter:: periodic(3) = [ .true., .false., .true. ]
      
      !> Control bounday condition
      integer,parameter:: periodicb=1,reflection=2,outflow=3
      integer,parameter:: boundary_xin=periodicb , boundary_xout=periodicb
      integer,parameter:: boundary_yin=reflection, boundary_yout=reflection
      integer,parameter:: boundary_zin=periodicb , boundary_zout=periodicb

      !> Control DATA-IO
      logical,parameter:: asciiout = .true. !! Ascii-files are additionaly damped.
      logical,parameter:: benchmarkmode = .true. !! If true, only initial and final outputs are damped. 

!$omp declare target (ncomp)      
      end module config
