
exed := ../exe

exe=Simulation.x

all: ${exed}/${exe}

.PHONY: all clean

#################
# C++
#################

cc := mpic++
copt := -std=c++17  -O3
#coptgpu := -mp=gpu -Minfo=mp,accel -gpu=cc80,mem:managed
coptgpu := -mp=gpu -Minfo=mp,accel -gpu=cc80
coptdebug := -g -traceback

# CPU environment for debug
#cc := icpx
#copt := -std=c++17  -O3 -lrt
#coptgpu := 
#coptdebug := -g -traceback

.SUFFIXES: .cpp .o

.cpp.o: 
	${cc} ${copt} ${coptgpu} ${coptdebug} -c $< 

OBJS = main.o mhd.o boundary.o output.o mpi_config.o

mpi_config.o : mpi_config.cpp mpi_config.hpp
output.o : output.cpp output.hpp
mhd.o: mhd.cpp mhd.hpp resolution.hpp mpi_config.o
boundary.o: boundary.cpp boundary.hpp mpi_config.o
main.o: main.cpp mhd.o boundary.o mpi_config.o output.o

${exed}/${exe}: ${OBJS}
	${cc} ${copt} ${coptgpu} ${coptdebug} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} ${ana} *.o *.mod *~
