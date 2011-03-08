IF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)
  SET(GoTools_FIND_QUIETLY TRUE)
ENDIF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)

FIND_PATH(GoTools_INCLUDE_DIRS
  NAMES GoTools/geometry/SplineSurface.h
  PATHS "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  /sima/libs/GoTools/include
)

FIND_LIBRARY(GoTools_LIBRARIES
  NAMES GoToolsCore
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  /sima/libs/GoTools/lib
  PATH_SUFFIXES GoTools
)

INCLUDE(FindPackageHandleStandardArgs)
IF(GoTools_LIBRARIES)
  find_package_handle_standard_args(GoTools DEFAULT_MSG
                                    GoTools_INCLUDE_DIRS GoTools_LIBRARIES)
ENDIF(GoTools_LIBRARIES)
