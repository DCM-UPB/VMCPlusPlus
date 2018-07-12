#!/bin/bash

source config.sh

\rm -f *.so
cd src/
\rm -f *.o *.so
echo "$CC $FLAGS $OPTFLAGS -fpic $IMCI $INFM -c *.cpp"
$CC $FLAGS $OPTFLAGS -fpic $IMCI $INFM -c *.cpp

case ${OS_NAME} in
    "Darwin")
        ROOT_FOLDER=$(dirname $(pwd))
        echo "$CC $FLAGS $OPTFLAGS -shared -install_name ${ROOT_FOLDER}/lib${LIBNAME}.so $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL"
        $CC $FLAGS $OPTFLAGS -shared -install_name ${ROOT_FOLDER}/lib${LIBNAME}.so $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL
        ;;
    "Linux")
        echo "$CC $FLAGS $OPTFLAGS -shared $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL"
        $CC $FLAGS $OPTFLAGS -shared $LMCI $LNFM $LGSL -o lib${LIBNAME}.so *.o $LIBMCI $LIBNFM $LIBGSL
        ;;
esac

mv lib${LIBNAME}.so ../
cd ..

echo
echo "Library ready!"
echo
echo "Help, how can I use it?"
echo "1)   $CC -I$(pwd)/src/ -c example.cpp"
echo "     $CC -L$(pwd) example.o -l${LIBNAME}"
echo "2)   $CC -I$(pwd)/src/ -L$(pwd) example.cpp -l${LIBNAME}"
