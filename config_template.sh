#!/bin/bash

OS_NAME=$(uname)

# Library name
LIBNAME="vmc"

# C++ compiler
CC="g++"

# MPI compiler wrapper
MPICC="mpic++"
USE_MPI=0

# C++ flags (std=c++11 is necessary)
FLAGS="-std=c++11 -Wall -Wno-unused-function"

# Optimization flags
OPTFLAGS="-O3 -flto"

# Debuggin flags
DEBUGFLAGS="-g -O0"

# GSL (GNU Scientific Library)
LGSL="-L/usr/local/lib"
LIBGSL="-lgsl -lgslcblas"

# MCIntegrator++ Library
MCI_FOLDER="/...../MCIntegratorPlusPlus"
IMCI="-I${MCI_FOLDER}/src/"
LMCI="-L${MCI_FOLDER}"
LIBMCI="-lmci"

#NoisyFunctionMinimization Library
NFM_FOLDER="/...../NoisyFunMin"
INFM="-I${NFM_FOLDER}/src/"
LNFM="-L${NFM_FOLDER}"
LIBNFM="-lnfm"
