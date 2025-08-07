
exed := ../exe

exe=Simulation.x

all: ${exed}/${exe}

.PHONY: all clean

#################
# Fortran
#################

fc= mpif90
fopt = -g -traceback -O2 -mcmodel=medium
foptgpu = -mp=gpu -Minfo=mp

.SUFFIXES: .f90 .o

.f90.o: 
	${fc} ${fopt} ${foptgpu} -c $< 

OBJS = main.o output.o boundary.o mhd.o mpiroutines.o

mpiroutines.o: mpiroutines.f90
mhd.o: mhd.f90 mpiroutines.o
boundary.o: boundary.f90 mhd.o mpiroutines.o 
output.o: output.f90 boundary.o mhd.o mpiroutines.o
main.o: main.f90 output.o boundary.o mhd.o mpiroutines.o

${exed}/${exe}: ${OBJS}
	${fc} ${fopt} ${foptgpu} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} *.o *.mod *~

