IF (SAMG_LIBRARIES)
  SET(SAMG_FIND_QUIETLY TRUE)
ENDIF(SAMG_LIBRARIES)

FIND_PATH(SAMG_INCLUDES
          NAMES
          samg.h
          PATHS
          $ENV{HOME}/include
)

FIND_LIBRARY(SAMG_LIBRARIES amg_coo 
             PATHS 
             $ENV{HOME}/lib
)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SAMG DEFAULT_MSG
                                  SAMG_LIBRARIES)

MARK_AS_ADVANCED(SAMG_LIBRARIES)
