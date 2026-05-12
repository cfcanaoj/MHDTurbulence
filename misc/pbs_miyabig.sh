#PBS -q regular-g
#PBS -l select=1:mpiprocs=1
#PBS -l walltime=00:10:00
#PBS -W group_list=xg25i092
#PBS -j oe

cd ${PBS_O_WORKDIR}
mpirun ./Simulation.x
