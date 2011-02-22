# - Tries to find the GoTools Trivariate library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoTrivariate_INCLUDE_DIRS ${GoTrivariate_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Trivariate header files")
  # Library
  SET(GoTrivariate_LIBRARIES GoTrivariate
    CACHE FILE "GoTools Trivariate library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoTrivariate_INCLUDE_DIRS 
  "GoTools/trivariate/SplineVolume.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "/usr/local/include"
  /sima/libs/GoTools/include
  "C:/Program Files (x86)/GoTools/include"
)

# Find library
FIND_LIBRARY(GoTrivariate_LIBRARIES
  NAMES GoTrivariate
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "/usr/local/lib"
  /sima/libs/GoTools/lib
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(GoTrivariate_FOUND FALSE)
IF(GoTrivariate_INCLUDE_DIRS AND GoTrivariate_LIBRARIES)
  SET(GoTrivariate_FOUND TRUE)
ENDIF(GoTrivariate_INCLUDE_DIRS AND GoTrivariate_LIBRARIES)
