exed := ../exe

exe = Simulation.x

all: ${exed}/${exe}

.PHONY: all clean preprocess

#################
# Fortran
#################

fc ?= mpifort
fopt ?= -g -traceback -O2

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SOLOMON_CPP ?= cc -E
else
  SOLOMON_CPP ?= cpp
endif

SOLOMON_DIR ?= ./solomon
USE_ACC_PARALLEL ?= 0
OPENACC_VERSION ?= 201711
OPENMP_VERSION ?= 202011
solomon_defs ?= -DOFFLOAD_BY_OPENACC

ifeq ($(USE_ACC_PARALLEL),1)
  solomon_defs += -DOFFLOAD_BY_OPENACC_PARALLEL
endif
foptoffload ?= -Minfo=accel -acc -mcmodel=medium

inc = -I$(SOLOMON_DIR)

SOURCES = main.F90 output.F90 boundary.F90 mhd.F90 mpiroutines.F90 config.F90
OBJS = $(SOURCES:.F90=.o)
GENS = $(addprefix spp/,$(SOURCES:.F90=.f90))

SOLOMON_HASH := \#
SOLOMON_FC = $(fc)
SOLOMON_FLAGS = $(fopt) $(foptoffload) $(solomon_defs) $(inc)
SOLOMON_PROBE_FILE = ._solomon.F
SOLOMON_PROBE_DEF := $(shell printf "$(SOLOMON_HASH)ifdef _OPENMP\n-D_OPENMP=_OPENMP\n$(SOLOMON_HASH)endif\n$(SOLOMON_HASH)ifdef _OPENACC\n-D_OPENACC=_OPENACC\n$(SOLOMON_HASH)endif\n" > $(SOLOMON_PROBE_FILE); $(SOLOMON_FC) $(SOLOMON_FLAGS) -E $(SOLOMON_PROBE_FILE) 2>/dev/null | grep -- "-D_"; rm -f $(SOLOMON_PROBE_FILE))
ifneq ($(findstring -DOFFLOAD_BY_OPENACC,$(solomon_defs)),)
  ifeq ($(findstring -D_OPENACC,$(SOLOMON_PROBE_DEF)),)
    SOLOMON_BACKEND_DEF += -D_OPENACC=$(OPENACC_VERSION)
  endif
endif
ifneq ($(findstring -DOFFLOAD_BY_OPENMP_TARGET,$(solomon_defs)),)
  ifeq ($(findstring -D_OPENMP,$(SOLOMON_PROBE_DEF)),)
    SOLOMON_BACKEND_DEF += -D_OPENMP=$(OPENMP_VERSION)
  endif
endif
SOLOMON_DEF = $(SOLOMON_PROBE_DEF) $(SOLOMON_BACKEND_DEF) $(filter -D% -I%,$(SOLOMON_FLAGS))

.SUFFIXES:
.SECONDARY: $(GENS)

preprocess: $(GENS)

spp/%.f90: %.F90
	mkdir -p spp
	$(SOLOMON_CPP) -P -x c++ -DSOLOMON_FORTRAN $(SOLOMON_DEF) $< > spp/$*.i.f90
	sed 's/^#pragma /!$$/g' spp/$*.i.f90 > $@
	rm -f spp/$*.i.f90

%.o: spp/%.f90
	${fc} ${fopt} ${foptoffload} -c $< -o $@

config.o: spp/config.f90
mpiroutines.o: spp/mpiroutines.f90 config.o
mhd.o: spp/mhd.f90 mpiroutines.o config.o
boundary.o: spp/boundary.f90 mhd.o mpiroutines.o config.o
output.o: spp/output.f90 boundary.o mhd.o mpiroutines.o
main.o: spp/main.f90 output.o boundary.o mhd.o mpiroutines.o

${exed}/${exe}: ${OBJS}
	mkdir -p ${exed}
	${fc} ${fopt} ${foptoffload} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} *.o *.mod *~
	rm -rf spp
