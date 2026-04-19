#! /bin/bash
#SBATCH --partition=dgx-full
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:4
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log

# usage sbatch sj_g00.sh
# other useful commands
# sinfo
# squeue

module load nvhpc/25.7
date >& out${SLURM_JOB_ID}.log
time mpiexec --bind-to none -n 4 ./Simulation.x >> out${SLURM_JOB_ID}.log
date >> out${SLURM_JOB_ID}.log

