exed := ../exe

exe = Simulation.x

all: ${exed}/${exe}

.PHONY: all clean preprocess

#################
# Fortran
#################

fc ?= mpifort
fopt ?= -g -traceback -O2
#foptopenacc ?= -Minfo=accel -acc -mcmodel=medium
foptopenacc ?= -Minfo=mp -mp=gpu -mcmodel=medium


cpp ?= cc -E -x c
#solomon_defs = -DSOLOMON_FORTRAN -DOFFLOAD_BY_OPENACC -D_OPENACC=201711 -I./solomon
solomon_defs = -DSOLOMON_FORTRAN -DOFFLOAD_BY_OPENMP_TARGET -DOFFLOAD_BY_OPENMP_TARGET_DISTRIBUTE -I./solomon

SOURCES = main.F90 output.F90 boundary.F90 mhd.F90 mpiroutines.F90 config.F90
OBJS = $(SOURCES:.F90=.o)
GENS = $(addprefix spp/,$(SOURCES:.F90=.f90))

.SUFFIXES:
.SECONDARY: $(GENS)

preprocess: $(GENS)

spp/%.f90: %.F90
	mkdir -p spp
	sed 's,//,__SOLOMON_FC_CONCAT__,g' $< > spp/$*.src.F90
	$(cpp) -P $(solomon_defs) spp/$*.src.F90 > spp/$*.i.f90
	sed -e 's/^#pragma /!$$/g' -e 's/__SOLOMON_FC_CONCAT__/\/\//g' spp/$*.i.f90 > $@
	rm -f spp/$*.src.F90 spp/$*.i.f90

%.o: spp/%.f90
	${fc} ${fopt} ${foptopenacc} -c $< -o $@

config.o: spp/config.f90
mpiroutines.o: spp/mpiroutines.f90 config.o
mhd.o: spp/mhd.f90 mpiroutines.o config.o
boundary.o: spp/boundary.f90 mhd.o mpiroutines.o config.o
output.o: spp/output.f90 boundary.o mhd.o mpiroutines.o
main.o: spp/main.f90 output.o boundary.o mhd.o mpiroutines.o

${exed}/${exe}: ${OBJS}
	mkdir -p ${exed}
	${fc} ${fopt} ${foptopenacc} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} *.o *.mod *~
	rm -rf spp
