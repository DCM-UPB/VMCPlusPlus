#!/bin/sh
# Scans the library using clang-tidy, with a large set of checks enabled
# First argument can be used to add extra arguments, like -fix
. ./config.sh
for file in src/*.*pp test/main.cpp test/*/*.*pp examples/*/*.*pp; do
    clang-tidy -p build $file $1 -checks=bugprone-*,cppcoreguidelines-*,clang-analyzer-*,google-*,llvm-*,misc-*,modernize-*,mpi-*,performance-*,portability-*,readability-*,-clang-analyzer-optin.performance.Padding,-cppcoreguidelines-owning-memory,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-cppcoreguidelines-pro-bounds-constant-array-index,-cppcoreguidelines-special-member-functions,-google-runtime-references,-google-build-using-namespace,-llvm-header-guard -header-filter=.* -- -Iinclude/ -Itest/common/ -I${MCI_ROOT}/include -I${NFM_ROOT}/include -DUSE_MPI=1
done
