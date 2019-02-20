#!/bin/sh
cd ../../build/examples
mpirun -np $1 ./ex7.exe
