
exed := ../exeacc

exe=Simulation.x

all: ${exed}/${exe}

.PHONY: all clean

#################
# Fortran
#################

fc= mpif90
fopt = -g -traceback -O2
foptopenacc = -Minfo=accel -acc -mcmodel=medium

#fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

.SUFFIXES: .f90 .o

.f90.o: 
	${fc} ${fopt} ${foptopenacc} -c $< 

OBJS = main.o output.o boundary.o mhd.o mpiroutines.o config.o

config.o: config.f90
mpiroutines.o: mpiroutines.f90 config.o
mhd.o: mhd.f90 mpiroutines.o config.o
boundary.o: boundary.f90 mhd.o mpiroutines.o config.o 
output.o: output.f90 boundary.o mhd.o mpiroutines.o
main.o: main.f90 output.o boundary.o mhd.o mpiroutines.o

${exed}/${exe}: ${OBJS}
	${fc} ${fopt} ${foptopenacc} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} *.o *.mod *~
