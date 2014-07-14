#!/bin/bash

export OMP_NUM_THREADS=8

mpiexec -np 128 ./vmc.out multiDir.def

