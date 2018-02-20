#!/bin/bash

# After using this script it is necessary to run again the build.sh script
# for generating again the library with the optimization flags

OS_NAME=$(uname)
source ../../../config.sh
DEBUGFLAGS="-g -O0"

\rm -f exe
\rm -f *.o

ROOT_FOLDER=$(dirname $(dirname $(dirname $(pwd))))


## Build the debugging main executable
$CC $FLAGS $DEBUGFLAGS -Wall $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp

case ${OS_NAME} in
    "Linux")
        $CC $FLAGS $FLAG_TO_USE $LMCI $LNFM $LFFNN $LGSL -I${ROOT_FOLDER}/src/ -L${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
        ;;
    "Darwin")
        $CC $FLAGS $FLAG_TO_USE -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
        ;;
    *)
        echo "The detected operating system is not between the known ones (Linux and Darwin)"
        ;;
esac

# Run the debugging executable
./exe
