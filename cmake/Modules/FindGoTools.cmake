# - Tries to find the GoTools Core library
#
# Written by: jan.b.thomassen@sintef.no
#


# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoTools_INCLUDE_DIRS ${GoToolsCore_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Core header files")
  # Library
  SET(GoTools_LIBRARIES GoToolsCore
    CACHE FILE "GoTools Core library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoTools_INCLUDE_DIRS "GoTools/geometry/SplineSurface.h"
  /sima/libs/GoTools/include
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  )


# Find library
FIND_LIBRARY(GoTools_LIBRARIES
  NAMES GoToolsCore 
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  /sima/libs/GoTools/lib
  PATH_SUFFIXES GoTools
  )


# Check that we have found everything
SET(GoTools_FOUND FALSE)
IF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)
  SET(GoTools_FOUND TRUE)
ENDIF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)
