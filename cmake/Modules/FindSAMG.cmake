find_path(
  SAMG_INCLUDES
  NAMES
    samg.h
  PATHS
    $ENV{HOME}/include
)

find_library(
  SAMG_LIBRARIES
  NAMES
    amg_coo
  PATHS
    $ENV{HOME}/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SAMG DEFAULT_MSG
                                  SAMG_INCLUDES SAMG_LIBRARIES)

mark_as_advanced(SAMG_INCLUDES SAMG_LIBRARIES)

if(SAMG_FOUND)
  add_library(SAMG::SAMG UNKNOWN IMPORTED)
  set_target_properties(SAMG::SAMG PROPERTIES IMPORTED_LOCATION ${SAMG_LIBRARIES})
  target_include_directories(SAMG::SAMG INTERFACE ${SAMG_INCLUDES})
endif()
