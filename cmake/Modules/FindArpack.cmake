
IF (Arpack_LIBRARIES)
  SET(Arpack_FIND_QUIETLY TRUE)
ENDIF (Arpack_LIBRARIES)

FIND_PATH(Arpack_INCLUDES
          NAMES
          arpack++/arpackf.h
          PATHS
          $ENV{HOME}/include
)

FIND_LIBRARY(Arpack_LIBRARIES arpack 
             PATHS 
             $ENV{HOME}/lib
)

IF(Arpack_LIBRARIES AND CMAKE_COMPILER_IS_GNUCXX)
  SET(Arpack_LIBRARIES ${Arpack_LIBRARIES})
ENDIF(Arpack_LIBRARIES AND CMAKE_COMPILER_IS_GNUCXX)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Arpack DEFAULT_MSG
                                  Arpack_LIBRARIES)

MARK_AS_ADVANCED(Arpack_LIBRARIES)
