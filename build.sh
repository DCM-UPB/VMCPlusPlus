#!/bin/sh

source ./config.sh
mkdir -p build
cd build
cmake -DUSER_CXX_FLAGS="${CXX_FLAGS}" -DUSE_MPI="${USE_MPI}" -DMCI_ROOT_DIR="${MCI_ROOT}" -DNFM_ROOT_DIR="${NFM_ROOT}" -DGSL_ROOT_DIR="${GSL_ROOT}" ..
make
cd ..
