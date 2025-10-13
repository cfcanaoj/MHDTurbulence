#!/bin/bash
#SBATCH --job-name=MHDturbulence
#SBATCH --partition=M-large-cfca
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -o out%j.log
#SBATCH -e err%j.log
#SBATCH --time=1:00:00

source /work/opt/local/bin/enable-oneapi.sh
module switch PrgEnv-oneapi/2024.2 PrgEnv-oneapi/2025.1

cd ${SLURM_SUBMIT_DIR}

# OpenMP
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY=compact

date >& out${SLURM_JOB_ID}.log
time srun -c ${SLURM_CPUS_PER_TASK} ./Simulation.x >> out${SLURM_JOB_ID}.log
date >> out${SLURM_JOB_ID}.log



