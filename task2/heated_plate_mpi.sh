#! /bin/bash
#

mpicc -O3 heated_plate_mpi.c -o heated_plate_mpi -lm
mpirun -np 4 heated_plate_mpi 
