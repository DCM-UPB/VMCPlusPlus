#!/bin/sh

. ./config.sh
mkdir -p build
cd build
cmake -DCMAKE_CXX_COMPILER="${CXX_COMPILER}" -DUSER_CXX_FLAGS="${CXX_FLAGS}" -DUSE_COVERAGE="${USE_COVERAGE}" -DUSE_MPI="${USE_MPI}" -DMCI_ROOT_DIR="${MCI_ROOT}" -DNFM_ROOT_DIR="${NFM_ROOT}" -DGSL_ROOT_DIR="${GSL_ROOT}" ..
make -j$(nproc)
