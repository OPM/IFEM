find_library(
  ARPACK_LIBRARIES
  NAMES
    arpack
  PATHS
    $ENV{HOME}/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG
                                  ARPACK_LIBRARIES)

if(ARPACK_FOUND)
  add_library(ARPACK::ARPACK UNKNOWN IMPORTED)
  set_target_properties(ARPACK::ARPACK PROPERTIES IMPORTED_LOCATION ${ARPACK_LIBRARIES})
endif()

mark_as_advanced(ARPACK_LIBRARIES)
