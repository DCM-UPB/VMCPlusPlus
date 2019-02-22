file(GLOB SOURCES "*.cpp")
add_library(vmc SHARED ${SOURCES})
target_link_libraries(vmc "${MCI_LIBRARY_DIR}" "${NFM_LIBRARY_DIR}" "${GSL_LIBRARIES}" "${MPI_CXX_LIBRARIES}") # shared libs
add_library(vmc_static STATIC ${SOURCES})
target_link_libraries(vmc_static "${MCI_STATIC_LIBRARY_DIR}" "${NFM_STATIC_LIBRARY_DIR}" "{$GSL_LIBRARIES}" "${MPI_CXX_LIBRARIES}") # static (+ some shared) libs
