#!/bin/bash

# use MPI for integration
USE_MPI=0

# C++ flags
CXX_FLAGS="-O3 -flto -Wall -Wno-unused-function"

# MCIntegrator++ Library
MCI_ROOT="/...../MCIntegratorPlusPlus"

#NoisyFunctionMinimization Library
NFM_ROOT="/...../NoisyFunMin"

# GNU Scientific Library
GSL_ROOT="" # provide a path if not in system location

# currently unused:
# C++ compiler
# CXX="g++"
