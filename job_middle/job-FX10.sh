#!/bin/bash
#PJM -N MVMC-FX10-FIPP-64
#PJM --rsc-list "node=64"
#PJM --rsc-list "elapse=0:30:00"
#PJM --mpi "proc=128"
#PJM --mpi "rank-map-bychip"
#PJM -j
#PJM -S

source /home/system/Env_base
set -x
date
hostname
cat /etc/system-release

WDIR=${HOME}/mvmc/mVMC-mini/job_middle
cd $WDIR; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
#	rm -rf prof_dir Lx10Ly10_J1.0/*.dat
export OMP_NUM_THREADS=8
fipp -C -d prof_dir -Icall,hwm -Puserfunc \
mpiexec -n 128 ./vmc.out multiDir.def

