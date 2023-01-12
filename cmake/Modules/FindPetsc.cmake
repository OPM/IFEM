find_package(PkgConfig)

set(OLD_PKG $ENV{PKG_CONFIG_PATH})
set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig)
set(OLD_ALLOW $ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS})
set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
pkg_check_modules(PETSC PETSc>=3.6.3)
set(ENV{PKG_CONFIG_PATH} ${OLD_PKG})
set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} ${OLD_ALLOW})

set(PETSC_LIBRARIES ${PETSC_STATIC_LDFLAGS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Petsc DEFAULT_MSG
                                  PETSC_INCLUDE_DIRS PETSC_LIBRARIES)
