#!/bin/bash

source config.sh

\rm -f *.so
cd src/
\rm -f *.o *.so
echo "$CC -D USE_MPI=$USE_MPI $FLAGS $OPTFLAGS -fpic $IMCI $INFM -c *.cpp"
$CC -D USE_MPI=$USE_MPI $FLAGS $OPTFLAGS -fpic $IMCI $INFM -c *.cpp

case ${OS_NAME} in
    "Darwin")
        ROOT_FOLDER=$(dirname $(pwd))
        echo "$CC -D USE_MPI=$USE_MPI $FLAGS $OPTFLAGS -shared -install_name ${ROOT_FOLDER}/lib${LIBNAME}.so $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL"
        $CC -D USE_MPI=$USE_MPI $FLAGS $OPTFLAGS -shared -install_name ${ROOT_FOLDER}/lib${LIBNAME}.so $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL
        ;;
    "Linux")
        echo "$CC -D USE_MPI=$USE_MPI $FLAGS $OPTFLAGS -shared $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL"
        $CC -D USE_MPI=$USE_MPI $FLAGS $OPTFLAGS -shared $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL
        ;;
esac

mv lib${LIBNAME}.so ../
cd ..

echo
echo "Library ready!"
echo
echo "Help, how can I use it?"

case ${USE_MPI} in
    "0")
        MYCC=$CC
        ;;
    "1")
        MYCC=$MPICC
        ;;
esac

echo "1)   $MYCC -I$(pwd)/src/ -c example.cpp"
echo "     $MYCC -L$(pwd) example.o -l${LIBNAME}"
echo "2)   $MYCC -I$(pwd)/src/ -L$(pwd) example.cpp -l${LIBNAME}"
