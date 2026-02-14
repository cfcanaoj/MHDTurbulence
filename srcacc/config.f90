      module config
      implicit none
      real(8),parameter:: timemax = 15.0d0
      real(8),parameter:: dtout = timemax/100

      integer,parameter:: nhymax = 600000
      integer,parameter:: nhydis = nhymax/100
      
      integer,parameter:: ngridtotal1 = 150
      integer,parameter:: ngridtotal2 = 150
      integer,parameter:: ngridtotal3 = 150

      real(8),parameter:: x1min = -0.5d0, x1max = 0.5d0
      real(8),parameter:: x2min = -1.0d0, x2max = 1.0d0
      real(8),parameter:: x3min = -0.5d0, x3max = 0.5d0      
      integer,parameter:: ncomp=1 ! kinds of composition
      
      integer,parameter:: ntiles(3) = [ 2,2,1 ]
      logical,parameter:: periodic(3) = [ .true., .false., .true. ]
      
      integer,parameter:: periodicb=1,reflection=2,outflow=3
      integer,parameter:: boundary_xin=periodicb , boundary_xout=periodicb
      integer,parameter:: boundary_yin=reflection, boundary_yout=reflection
      integer,parameter:: boundary_zin=periodicb , boundary_zout=periodicb
      
      logical,parameter:: asciiout= .true.

      end module config
