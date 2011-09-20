IF(LRSpline_INCLUDE_DIRS AND LRSpline_LIBRARIES)
  SET(LRSpline_FIND_QUIETLY TRUE)
ENDIF(LRSpline_INCLUDE_DIRS AND LRSpline_LIBRARIES)

UNSET(LRSpline_INCLUDE_DIRS CACHE)
UNSET(LRSpline_LIBRARIES CACHE)

FIND_PATH(LRSpline_INCLUDE_DIRS
  NAMES LRSpline/LRSplineSurface.h
  PATHS "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  /sima/libs/LRSpline/include
)

FIND_LIBRARY(LRSpline_LIBRARIES
  NAMES LRSpline
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  /sima/libs/LRSpline/lib
  PATH_SUFFIXES LRSpline
)

INCLUDE(FindPackageHandleStandardArgs)
IF(LRSpline_LIBRARIES)
  find_package_handle_standard_args(LRSpline DEFAULT_MSG
                                    LRSpline_INCLUDE_DIRS LRSpline_LIBRARIES)
ENDIF(LRSpline_LIBRARIES)

