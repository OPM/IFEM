IF (ARPACK_LIBRARIES)
  SET(ARPACK_FIND_QUIETLY TRUE)
ENDIF (ARPACK_LIBRARIES)

FIND_LIBRARY(ARPACK_LIBRARIES
  NAMES arpack
  PATHS $ENV{HOME}/lib
  /sima/libs/ARPACK
 # For kongull until they get their act together
  /share/apps/modulessoftware/arpack/96/lib
)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG
                                  ARPACK_LIBRARIES)

MARK_AS_ADVANCED(ARPACK_LIBRARIES)
