#!/bin/bash

export OMP_NUM_THREADS=1

mpiexec -np 1 ./vmc.out multiDir.def

