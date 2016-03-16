#!/bin/bash

# After using this script it is necessary to run again the build.sh script
# for generating again the library with the optimization flags

source ../config.sh
DEBUGFLAGS="-g -O0"

\rm -f exe
\rm -f *.o
\rm -f ../src/*.o
\rm -f ../*.so

#runtime dynamic library path
RPATH="$(pwd)/..:${MCI_FOLDER}:${NFM_FOLDER}"

# Build the library using the debugging flags
cd ../src
   $CC $DEBUGFLAGS -std=c++11 $IMCI $INFM -Wall -Werror -fpic -c *.cpp
   $CC $DEBUGFLAGS -std=c++11 $LMCI $LNFM -shared -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM
   mv lib*.so ../
cd ../debug
echo "Rebuilt the library with the debugging flags"

# Build the debugging main executable
echo "$CC $FLAGS $DEBUGFLAGS $IMCI $INFM -I$(pwd)/../src/ -c *.cpp"
$CC $FLAGS $DEBUGFLAGS -Wall $IMCI $INFM -I$(pwd)/../src/ -c *.cpp
echo "$CC $FLAGS $DEBUGFLAGS $LMCI $LNFM -I$(pwd)/../src -L$(pwd)/.. -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}" $LIBMCI $LIBNFM
$CC $FLAGS $DEBUGFLAGS $LMCI $LNFM -I$(pwd)/../src/ -L$(pwd)/../ -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} $LIBMCI $LIBNFM
echo "Rebuilt the debugging executable"
echo ""
echo ""

# Run the debugging executable
valgrind --leak-check=full --track-origins=yes ./exe


