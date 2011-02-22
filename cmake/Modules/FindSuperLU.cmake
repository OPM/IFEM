IF(SuperLU_INCLUDES AND SuperLU_LIBRARIES)
  SET(SUPERLU_FIND_QUIETLY TRUE)
ENDIF(SuperLU_INCLUDES AND SuperLU_LIBRARIES)

FIND_PATH(SuperLU_INCLUDES
          NAMES
          supermatrix.h
          PATHS
          $ENV{HOME}/include
          /sima/libs/SuperLU_4.0/include
          PATH_SUFFIXES superlu
)

FIND_PATH(SuperLU_MT_INCLUDES
          NAMES
          slu_mt_machines.h
          PATHS
          $ENV{HOME}/include
          /sima/libs/SuperLU_4.0/include
          PATH_SUFFIXES superlu
)

FIND_LIBRARY(SuperLU_LIBRARIES superlu 
             PATHS 
             /sima/libs/SuperLU_4.0/lib
             $ENV{HOME}/lib
)

FIND_LIBRARY(SuperLU_MT_LIBRARIES superlu_mt
             PATHS 
             /sima/libs/SuperLU_4.0/lib
             $ENV{HOME}/lib
)

IF(SuperLU_LIBRARIES AND CMAKE_COMPILER_IS_GNUCXX)
  SET(SuperLU_LIBRARIES ${SuperLU_LIBRARIES} -lgfortran)
ENDIF(SuperLU_LIBRARIES AND CMAKE_COMPILER_IS_GNUCXX)

IF(SuperLU_MT_LIBRARIES AND CMAKE_COMPILER_IS_GNUCXX)
  SET(SuperLU_MT_LIBRARIES ${SuperLU_MT_LIBRARIES} -lgfortran)
ENDIF(SuperLU_MT_LIBRARIES AND CMAKE_COMPILER_IS_GNUCXX)

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
