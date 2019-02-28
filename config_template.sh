#!/bin/sh

#C++ compiler
CXX_COMPILER="g++"

# C++ flags
CXX_FLAGS="-O3 -flto -march=native -Wall -Wno-unused-function"

# add coverage flags
USE_COVERAGE=0

# use MPI for integration
USE_MPI=0

# MCIntegrator++ Library
MCI_ROOT="/...../MCIntegratorPlusPlus"

# NoisyFunctionMinimization Library
NFM_ROOT="/...../NoisyFunMin"

# GNU Scientific Library
GSL_ROOT="" # provide a path if not in system location
