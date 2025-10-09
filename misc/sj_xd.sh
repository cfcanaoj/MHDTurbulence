#!/bin/bash
#SBATCH --job-name=MHDturbulence
#SBATCH --partition=M-large-cfca
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -o slmlog%J.out
#SBATCH --time=1:00:00

source /work/opt/local/bin/enable-oneapi.sh
module switch PrgEnv-oneapi/2024.2 PrgEnv-oneapi/2025.1

cd ${SLURM_SUBMIT_DIR}

# OpenMP
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY=compact

date >& log.$SLURM_JOB_ID
time srun -c ${SLURM_CPUS_PER_TASK} ./Simulation.x >> log.$SLURM_JOB_ID
date >> log.$SLURM_JOB_ID



