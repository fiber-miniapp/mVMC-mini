#!/bin/bash
#PJM -N MVMC-K-128
#PJM --rsc-list "node=128"
#PJM --rsc-list "elapse=2:00:00"
#PJM --mpi "proc=128"
#PJM --mpi "rank-map-bychip"
# stage io files
#PJM --stg-transfiles all
#PJM --stgin-basedir "/home/ra000004/a03155/mvmc/mVMC-mini"
#PJM --stgout-basedir "/home/ra000004/a03155/mvmc/scripts/logs_K/middle-128"
# rank directories
#PJM --mpi "use-rankdir"
#PJM --stgin "rank=* src/vmc.out %r:vmc.out"
#PJM --stgin "rank=* job_middle/multiDir.def %r:multiDir.def"
#PJM --stgin-dir "rank=* job_middle/Lx10Ly10_J1.0 %r:Lx10Ly10_J1.0"
#PJM --stgout-dir "rank=* %r:Lx10Ly10_J1.%r ./Lx10Ly10_J1.%r"
#PJM -j
#PJM -S

source /work/system/Env_base
set -x
date
hostname
cat /etc/system-release
pwd

export OMP_NUM_THREADS=8
time -p mpiexec -n 128 ./vmc.out multiDir.def

sleep 60
ls -goR ../
cat Lx10Ly10_J1.0/zvo_out_000.dat

