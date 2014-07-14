#!/bin/bash
#BSUB -J MVMC-M128-V0
#BSUB -o MVMC-M128-V0-%J
#BSUB -n 128
#BSUB -R "span[ptile=2]"
#BSUB -x
module load intel impi mkl
module list
set -x
date
hostname
TEST_DIR=${HOME}/mvmc/mVMC-mini/job_middle
cd ${TEST_DIR}
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi

export OMP_NUM_THREADS=8

mpirun -np 128 ./vmc.out multiDir.def

