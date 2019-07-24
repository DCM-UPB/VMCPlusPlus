#!/bin/sh
cd ../../build/examples
if [[ "$1" == "" ]]; then
    echo "!----------------------------------------------------------------------!"
    echo "! Warning: Number of threads not specified -> running single-threaded. !"
    echo "! To use parallel MC integration, pass the #threads: e.g. ./run.sh 2   !"
    echo "!----------------------------------------------------------------------!"
    echo
    ./ex_mpi.exe
else
    mpirun --oversubscribe -np $1 ./ex_mpi.exe
fi;
