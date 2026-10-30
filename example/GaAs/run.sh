#!/bin/bash
#SBATCH --job-name=ex_GaAs
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=00:30:00
#SBATCH --output=eo.%j
#SBATCH --error=eo.%j
# dev QE example, GaAs: plain PBE ground state with dev pw.x, then dev hp.x
# one-shot ACBN0 (acbn0 = 1, on-site U) and eACBN0 (acbn0 = 2, U + V).
source /etc/profile.d/lmod.sh
module purge
module load intel/19.1.3.304
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AUTO
BIN=/home/yjkang/projects/QE/qe-7.3_thiswork/bin
PWX=$BIN/pw.x.e2801ae      # PW/src identical to official QE 7.3
HPX=$BIN/hp.x.a1ad49f      # dev hp.x, commit a1ad49f
NPOOL=32
cd $SLURM_SUBMIT_DIR
md5sum $PWX $HPX
for f in scf.out hp_U.out hp_UV.out out; do [ -e $f ] && { echo "$f exists, refusing to overwrite"; exit 2; }; done
mkdir -p out
srun --mpi=pmi2 --kill-on-bad-exit=1 $PWX -npool $NPOOL -inp scf.in   > scf.out   2>&1 ; echo "pw.x exit $?"
srun --mpi=pmi2 --kill-on-bad-exit=1 $HPX -npool $NPOOL -inp hp_U.in  > hp_U.out  2>&1 ; echo "hp.x U exit $?"
mv out/eACBN0_hub.dat eACBN0_hub_U.dat
srun --mpi=pmi2 --kill-on-bad-exit=1 $HPX -npool $NPOOL -inp hp_UV.in > hp_UV.out 2>&1 ; echo "hp.x U+V exit $?"
mv out/eACBN0_hub.dat eACBN0_hub_UV.dat
