find_package(PkgConfig)

set(OLD_PKG $ENV{PKG_CONFIG_PATH})
set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig)
pkg_check_modules(ISTL dune-istl)
pkg_check_modules(ISTLc dune-common)
list(APPEND ISTL_LIBRARIES ${ISTLc_LIBRARIES})
set(ENV{PKG_CONFIG_PATH} ${OLD_PKG})

add_library(Dune::ISTL INTERFACE IMPORTED)
target_link_libraries(Dune::ISTL INTERFACE ${ISTL_LIBRARIES})
target_compile_definitions(Dune::ISTL INTERFACE HAVE_NULLPTR=1
                                                ISTL_VERSION=\"${ISTL_VERSION}\"
                                                HAS_ISTL=1)
string(REPLACE "." ";" VERSION_LIST "${ISTL_VERSION}")
if(VERSION_LIST)
  list(GET VERSION_LIST 0 DUNE_ISTL_VERSION_MAJOR)
  list(GET VERSION_LIST 1 DUNE_ISTL_VERSION_MINOR)
  list(LENGTH VERSION_LIST nver)
  if(nver GREATER 2)
    list(GET VERSION_LIST 2 DUNE_ISTL_VERSION_PATCH)
  else()
    set(DUNE_ISTL_VERSION_PATCH 0)
  endif()
  target_compile_definitions(
    Dune::ISTL
    INTERFACE
      DDUNE_ISTL_VERSION_MAJOR=${DUNE_ISTL_VERSION_MAJOR}
      DUNE_ISTL_VERSION_MINOR=${DUNE_ISTL_VERSION_MINOR}
      DUNE_ISTL_VERSION_REVISION=${DUNE_ISTL_VERSION_PATCH}
  )
endif()

if(TARGET SuperLU::SuperLU)
  target_compile_definitions(
    Dune::ISTL
    INTERFACE
      HAVE_SUPERLU=1
      SUPERLU_POST_2005_VERSION=1
      SUPERLU_NTYPE=1
      SUPERLU_MIN_VERSION_4_3=1
      SUPERLU_INT_TYPE=int
      SUPERLU_MIN_VERSION_5=1
  )
  target_link_libraries(Dune::ISTL INTERFACE SuperLU::SuperLU)
endif()

if(TARGET SuiteSparse::umfpack)
  target_link_libraries(Dune::ISTL INTERFACE SuiteSparse::umfpack)
  target_compile_definitions(
    Dune::ISTL
    INTERFACE
      HAVE_UMFPACK=1
      HAVE_SUITESPARSE_UMFPACK=1
  )
endif()

target_include_directories(Dune::ISTL INTERFACE ${ISTL_INCLUDEDIR})
list(APPEND ISTL_INCLUDE_DIRS ${ISTL_INCLUDEDIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ISTL DEFAULT_MSG
                                  ISTL_INCLUDE_DIRS)

mark_as_advanced(ISTL_INCLUDE_DIRS ISTL_INCLUDEDIR ISTL_LIBRARIES)
