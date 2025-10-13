#! /bin/bash
#SBATCH --partition=dgx-full
#SBATCH --nodes=1
#SBATCH --gpu-bind=closest
#SBATCH --ntasks=4
#SBATCH --gres=gpu:4
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log

# usage sbatch sj_g00.sh
# other useful commands
# sinfo
# squeue

module load nvhpc
date >& out${SLURM_JOB_ID}.log
time mpiexec -n 4 ./Simulation.x >> out${SLURM_JOB_ID}.log
date >> out${SLURM_JOB_ID}.log

