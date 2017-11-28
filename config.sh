#!/bin/bash

# Library name
LIBNAME="vmc"

# C++ compiler
CC="g++"

# C++ flags (std=c++11 is necessary)
FLAGS="-std=c++11 -Wall"

# Optimization flags
OPTFLAGS="-O3"

# Debuggin flags
DEBUGFLAGS="-g -O0"

# MCIntegrator++ Library
MCI_FOLDER="/Users/francesco086/Code/MCIntegrator++"
IMCI="-I${MCI_FOLDER}/src/"
LMCI="-L${MCI_FOLDER}"
LIBMCI="-lmci"

#NoisyFunctionMinimization Library
NFM_FOLDER="/Users/francesco086/Code/NoisyFunOpt"
INFM="-I${NFM_FOLDER}/src/"
LNFM="-L${NFM_FOLDER}"
LIBNFM="-lnfm"

#FFNN Library (used in ex3)
FFNN_FOLDER="/Users/francesco086/Code/FFNN"
INFM="-IFFNN_FOLDER/src/"
LNFM="-LFFNN_FOLDER"
LIBNAMEFFNN="ffnn"
LIBNFM="-lffnn"