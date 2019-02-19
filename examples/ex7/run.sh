#!/bin/bash
cd ../../build/examples
mpirun -np $1 ./ex7.exe
