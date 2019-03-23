#!/bin/sh
cd ../../build/examples
mpirun --oversubscribe -np $1 ./ex7.exe
