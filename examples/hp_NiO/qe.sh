#!/bin/bash
#SBATCH --job-name=NiO-hp
#SBATCH --partition=normal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --time=10:00:00
#SBATCH --output=eo.%j
#SBATCH --error=eo.%j


# qe.sh is executed by Slurm (not sourced), so `source` here is a normal,
# safe call into env.sh - its own internal guard only trips if env.sh
# itself gets run directly, which this isn't doing.
source $HOME/.config/qe/env.sh

export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AUTO

# total MPI ranks = nodes * ntasks-per-node = 2 * 64 = 128
# scf.in's 8x8x8 mesh reduces to 65 irreducible k-points (checked via a
# probe run, not the raw 8*8*8=512) - NPOOL=128 would exceed that and
# pw.x aborts with "some pools have no k-points". 64 is the largest
# divisor of 128 that's <= 65: 63 pools get 1 k-point, 1 pool gets 2,
# and each pool still gets 2 ranks for plane-wave/band parallelism.
# Retune if you change the k-point mesh.
NPOOL=64

srun --mpi=pmi2 qe_aCBN0-7.3_pw.x -npool $NPOOL -inp scf.in > scf.out
srun --mpi=pmi2 hp.x -npool $NPOOL -inp scf.in > scf.out

