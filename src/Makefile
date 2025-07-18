
exed := ../exe

exe=Simulation.x
ana=Analysis.x

all: ${exed}/${exe}

.PHONY: all clean

#################
# Fortran
#################

fc= mpif90
foptopenacc = -Minfo=accel -acc -mcmodel=medium
fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

.SUFFIXES: .f90 .o

.f90.o: 
	${fc} ${fopt} ${foptopenacc} -c $< 

OBJS = main.o mpiroutines.o

mpiroutines.o:
main.o: main.f90  mpiroutines.o 

${exed}/${exe}: ${OBJS}
	${fc} ${fopt} ${foptopenacc} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} ${ana} *.o *.mod *~
