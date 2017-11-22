#!/bin/bash

# After using this script it is necessary to run again the build.sh script
# for generating again the library with the optimization flags

source ../config.sh
OS_NAME=$(uname)
echo "System detected: $OS_NAME"

\rm -f exe
\rm -f *.o
\rm -f ../src/*.o
\rm -f ../*.so

# Build the library using the debugging flags
echo "Build the library using the debugging flags . . ."
cd ..
   LIBFOLDER=$(pwd)
   ./build_debug_library.sh
cd debug
echo "Rebuilt the library with the debugging flags"
echo

#runtime dynamic library path
RPATH="${LIBFOLDER}:${MCI_FOLDER}:${NFM_FOLDER}"

# Build the debugging main executable
echo "$CC $FLAGS $DEBUGFLAGS $IMCI $INFM -I$LIBFOLDER/src/ -c *.cpp"
$CC $FLAGS $DEBUGFLAGS -Wall $IMCI $INFM -I$LIBFOLDER/src/ -c *.cpp

# For Mac OS, the install name is wrong and must be corrected
case ${OS_NAME} in
   "Darwin")
      echo "$CC $FLAGS $DEBUGFLAGS -L$LIBFOLDER $LMCI $LNFM -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM"
      $CC $FLAGS $DEBUGFLAGS -L$LIBFOLDER $LMCI $LNFM -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM
      
      echo "install_name_tool -change lib${LIBNAME}.so $(pwd)/../lib${LIBNAME}.so exe"
      install_name_tool -change lib${LIBNAME}.so $(pwd)/../lib${LIBNAME}.so exe
      ;;
   "Linux")
      echo "$CC $FLAGS $DEBUGFLAGS $LMCI $LNFM -I$(pwd)/../src -L$(pwd)/.. -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}" $LIBMCI $LIBNFM
      $CC $FLAGS $DEBUGFLAGS $LMCI $LNFM -I$(pwd)/../src/ -L$(pwd)/../ -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} $LIBMCI $LIBNFM
      ;;
esac

echo "Rebuilt the debugging executable"
echo ""
echo ""




# Run the debugging executable
valgrind --leak-check=full --track-origins=yes ./exe
# ./exe


