IF(SuperLU_INCLUDES AND SuperLU_LIBRARIES)
  SET(SuperLU_FIND_QUIETLY TRUE)
ENDIF(SuperLU_INCLUDES AND SuperLU_LIBRARIES)

FIND_PATH(SuperLU_INCLUDES
  NAMES slu_ddefs.h
  PATHS $ENV{HOME}/include
  /sima/libs/SuperLU/include
  PATH_SUFFIXES superlu
)

FIND_PATH(SuperLU_MT_INCLUDES
  NAMES pdsp_defs.h
  PATHS $ENV{HOME}/include
  /sima/libs/SuperLU_MT/include
  PATH_SUFFIXES superlu_mt superlu
)

FIND_LIBRARY(SuperLU_LIBRARIES
  NAMES superlu
  PATHS $ENV{HOME}/lib
  /sima/libs/SuperLU/lib
)

FIND_LIBRARY(SuperLU_MT_LIBRARIES
  NAMES superlu_mt
  PATHS $ENV{HOME}/lib
  /sima/libs/SuperLU_MT/lib
)

INCLUDE(FindPackageHandleStandardArgs)
IF(SuperLU_LIBRARIES)
  find_package_handle_standard_args(SuperLU DEFAULT_MSG
                                    SuperLU_INCLUDES SuperLU_LIBRARIES)
ENDIF(SuperLU_LIBRARIES)
IF(SuperLU_MT_LIBRARIES)
  find_package_handle_standard_args(SuperLU DEFAULT_MSG
                                    SuperLU_MT_INCLUDES SuperLU_MT_LIBRARIES)
ENDIF(SuperLU_MT_LIBRARIES)

MARK_AS_ADVANCED(SuperLU_INCLUDES SuperLU_LIBRARIES
                 SuperLU_MT_INCLUDES SuperLU_MT_LIBRARIES)
