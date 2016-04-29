find_package(PkgConfig)

set(OLD_PKG $ENV{PKG_CONFIG_PATH})
set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig)
pkg_check_modules(ISTL dune-istl)
set(ENV{PKG_CONFIG_PATH} ${OLD_PKG})

list(APPEND ISTL_DEFINITIONS -DHAVE_NULLPTR=1 -DISTL_VERSION="${ISTL_VERSION}")
list(APPEND ISTL_INCLUDE_DIRS ${ISTL_INCLUDEDIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ISTL DEFAULT_MSG
                                  ISTL_INCLUDE_DIRS ISTL_LIBRARIES)
