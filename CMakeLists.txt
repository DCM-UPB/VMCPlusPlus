cmake_minimum_required(VERSION 3.5)
include(FindPackageHandleStandardArgs)

project(vmc LANGUAGES CXX VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USER_CXX_FLAGS}")

if (USE_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif ()

# find packages
message(STATUS "Configured MCI_ROOT_DIR: ${MCI_ROOT_DIR}")
message(STATUS "Configured NFM_ROOT_DIR: ${NFM_ROOT_DIR}")
message(STATUS "Configured GSL_ROOT_DIR: ${GSL_ROOT_DIR}")

if (USE_MPI)
    find_package(MPI)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI=1")
    message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
    message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
endif ()

find_package(GSL)
message(STATUS "GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIRS}")
message(STATUS "GSL_LIBRARIES: ${GSL_LIBRARIES}")

find_path(MCI_INCLUDE_DIR mci/MCIntegrator.hpp HINTS "${MCI_ROOT_DIR}/include/")
find_library(MCI_LIBRARY_DIR mci HINTS "${MCI_ROOT_DIR}/build/src/")
find_library(MCI_STATIC_LIBRARY_DIR mci_static HINTS "${MCI_ROOT_DIR}/build/src/")
find_package_handle_standard_args(libmci DEFAULT_MSG MCI_LIBRARY_DIR MCI_INCLUDE_DIR)
find_package_handle_standard_args(libmci_static DEFAULT_MSG MCI_STATIC_LIBRARY_DIR MCI_INCLUDE_DIR)

find_path(NFM_INCLUDE_DIR nfm/NoisyFunMin.hpp HINTS "${NFM_ROOT_DIR}/include")
find_library(NFM_LIBRARY_DIR nfm HINTS "${NFM_ROOT_DIR}/build/src/")
find_library(NFM_STATIC_LIBRARY_DIR nfm_static HINTS "${NFM_ROOT_DIR}/build/src/")
find_package_handle_standard_args(libnfm DEFAULT_MSG NFM_LIBRARY_DIR NFM_INCLUDE_DIR)
find_package_handle_standard_args(libnfm_static DEFAULT_MSG NFM_STATIC_LIBRARY_DIR NFM_INCLUDE_DIR)


message(STATUS "Configured CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "Configured CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# set header / library paths
include_directories(include/ "${MCI_INCLUDE_DIR}" "${NFM_INCLUDE_DIR}" "${GSL_INCLUDE_DIRS}" ${MPI_INCLUDE_PATH}) # headers

enable_testing()

# continue with subdirectories
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)
