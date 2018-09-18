#!/bin/bash

source ../../config.sh

OS_NAME=$(uname)

FLAG_TO_USE="${OPTFLAGS}"

case ${USE_MPI} in
    "0")
        MYCC=$CC
        ;;
    "1")
        MYCC=$MPICC
        ;;
esac

\rm -f exe
\rm -f *.o

# project root directory
ROOT_FOLDER=$(dirname $(dirname $(pwd)))

#runtime dynamic library path
RPATH="${ROOT_FOLDER}:${MCI_FOLDER}:${NFM_FOLDER}:${FFNN_FOLDER}"

# Build the debugging main executable
echo "$MYCC $FLAGS $FLAG_TO_USE $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp"
$MYCC $FLAGS $FLAG_TO_USE -Wall $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp

# For Mac OS, the install name is wrong and must be corrected
case ${OS_NAME} in
    "Darwin")
        echo "$MYCC $FLAGS $FLAG_TO_USE -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBFFNN $LIBGSL"
        $MYCC $FLAGS $FLAG_TO_USE -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
        ;;
    "Linux")
        echo "$MYCC $FLAGS $FLAG_TO_USE $LMCI $LNFM $LFFNN $LGSL -I${ROOT_FOLDER}/src -L${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}" $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
        $MYCC $FLAGS $FLAG_TO_USE $LMCI $LNFM $LFFNN $LGSL -I${ROOT_FOLDER}/src/ -L${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
        ;;
esac

echo "Rebuilt the executable file"
echo ""
echo ""

# Run the debugging executable
echo "Ready to run!"
echo ""
echo "--------------------------------------------------------------------------"
echo ""
echo ""
echo ""
./exe
# valgrind --leak-check=full --track-origins=yes ./exe
