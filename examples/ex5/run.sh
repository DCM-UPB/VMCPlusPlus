#!/bin/bash

source ../../config.sh

OS_NAME=$(uname)

\rm -f exe
\rm -f *.o

# project root directory
ROOT_FOLDER=$(dirname $(dirname $(pwd)))

#runtime dynamic library path
RPATH="${ROOT_FOLDER}:${MCI_FOLDER}:${NFM_FOLDER}"

# Build the debugging main executable
echo "$CC $FLAGS -headerpad_max_install_names $OPTFLAGS $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp"
$CC $FLAGS -headerpad_max_install_names $OPTFLAGS -Wall $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp

# For Mac OS, the install name is wrong and must be corrected
case ${OS_NAME} in
   "Darwin")
      echo "$CC $FLAGS -headerpad_max_install_names $OPTFLAGS -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBGSL $LIBFFNN"
      $CC $FLAGS -headerpad_max_install_names $OPTFLAGS -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
      
      echo "install_name_tool -change lib${LIBNAME}.so ${ROOT_FOLDER}/lib${LIBNAME}.so exe"
      install_name_tool -change lib${LIBNAME}.so ${ROOT_FOLDER}/lib${LIBNAME}.so exe
      
      echo "install_name_tool -change lib${LIBNAMEFFNN}.so ${FFNN_FOLDER}/lib${LIBNAMEFFNN}.so exe"
      install_name_tool -change lib${LIBNAMEFFNN}.so ${FFNN_FOLDER}/lib${LIBNAMEFFNN}.so exe
      ;;
   "Linux")
      echo "$CC $FLAGS $OPTFLAGS $LMCI $LNFM -I${ROOT_FOLDER}/src -L$${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}" $LIBMCI $LIBNFM
      $CC $FLAGS $OPTFLAGS $LMCI $LNFM -I${ROOT_FOLDER}/src/ -L${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} $LIBMCI $LIBNFM
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
