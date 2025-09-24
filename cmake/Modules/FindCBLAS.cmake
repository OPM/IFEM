find_package(CBLAS CONFIG QUIET)

if(NOT CBLAS_FOUND)
  find_package(BLAS)
  include(CheckFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
  check_function_exists(cblas_dgemm CBLAS_FOUND)
  if (NOT CBLAS_FOUND)
    find_library(
      CBLAS_LIBRARIES
      NAMES
       cblas
      PATHS
        $ENV{HOME}/lib
        /usr/lib/atlas
    )
  endif()
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(CBLAS DEFAULT_MSG
                                    CBLAS_LIBRARIES)
  if(CBLAS_FOUND)
    add_library(CBLAS::CBLAS UNKNOWN IMPORTED)
    set_target_properties(CBLAS::CBLAS PROPERTIES IMPORTED_LOCATION ${CBLAS_LIBRARIES})
  endif()
endif()

target_compile_definitions(CBLAS::CBLAS INTERFACE USE_CBLAS=1)

mark_as_advanced(CBLAS_LIBRARIES)
