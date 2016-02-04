IF(SuperLU_INCLUDES AND SuperLU_LIBRARIES)
  SET(SuperLU_FIND_QUIETLY TRUE)
ENDIF(SuperLU_INCLUDES AND SuperLU_LIBRARIES)

FIND_PATH(SuperLU_INCLUDES
  NAMES slu_ddefs.h
  PATHS $ENV{HOME}/include
  /sima/libs/SuperLU/include
  PATH_SUFFIXES superlu
)
FIND_LIBRARY(SuperLU_LIBRARIES
  NAMES superlu
  PATHS $ENV{HOME}/lib
  /sima/libs/SuperLU/lib
)

FIND_PATH(SuperLU_MT_INCLUDES
  NAMES slu_mt_ddefs.h
  PATHS $ENV{HOME}/include
  /sima/libs/SuperLU_MT/include
  PATH_SUFFIXES superlu_mt superlu
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
  find_package_handle_standard_args(SuperLU_MT DEFAULT_MSG
                                    SuperLU_MT_INCLUDES SuperLU_MT_LIBRARIES)
ENDIF(SuperLU_MT_LIBRARIES)

include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES ${SuperLU_LIBRARIES} ${CBLAS_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${SuperLU_INCLUDES} ${CBLAS_INCLUDES})
check_cxx_source_compiles("
    #include <slu_ddefs.h>

      int main()
      {
        dgssvx(NULL, NULL, NULL, NULL, NULL,
               NULL, NULL, NULL, NULL, NULL,
               NULL, 1, NULL, NULL, NULL,
               NULL, NULL, NULL, NULL,
               NULL, NULL, NULL);
      }
" HAVE_SUPERLU_5)

if(HAVE_SUPERLU_5)
  set(SuperLU_DEFINITIONS -DSUPERLU_VERSION=5)
else()
  set(SuperLU_DEFINITIONS -DSUPERLU_VERSION=4)
endif()

MARK_AS_ADVANCED(SuperLU_INCLUDES SuperLU_LIBRARIES SUPERLU_DEFINITIONS
                 SuperLU_MT_INCLUDES SuperLU_MT_LIBRARIES)
