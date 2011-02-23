IF(IFEM_INCLUDES AND IFEM_LIBRARIES)
  SET(IFEM_FIND_QUIETLY TRUE)
ENDIF(IFEM_INCLUDES AND IFEM_LIBRARIES)

IF(NOT DEFINED FORCE_SYSTEM_IFEM OR NOT ${FORCE_SYSTEM_IFEM})
  FIND_PATH(IFEM_INCLUDES
            NAMES
            SIMbase.h
            PATHS
            ${PROJECT_SOURCE_DIR}/../../src/
            PATH_SUFFIXES SIM)
ENDIF(NOT DEFINED FORCE_SYSTEM_IFEM OR NOT ${FORCE_SYSTEM_IFEM})

# Build is in-tree
IF(IFEM_INCLUDES)
  SET(IFEM_INCLUDES
      ${PROJECT_SOURCE_DIR}/../../src/ASM
      ${PROJECT_SOURCE_DIR}/../../src/Eig
      ${PROJECT_SOURCE_DIR}/../../src/Integrands
      ${PROJECT_SOURCE_DIR}/../../src/LinAlg
      ${PROJECT_SOURCE_DIR}/../../src/SIM
      ${PROJECT_SOURCE_DIR}/../../src/Utility)
  MESSAGE(STATUS "Using in-tree libIFEM")
  FIND_LIBRARY(IFEM_LIBRARIES 
               IFEM
               PATHS
               ${PROJECT_SOURCE_DIR}/../../${CMAKE_BUILD_TYPE}/lib
               ${PROJECT_SOURCE_DIR}/../../lib
               NO_DEFAULT_PATH)
  IF(NOT IFEM_LIBRARIES)
    MESSAGE(WARNING "Could not find the in-tree libIFEM library, we assume "
                    "it will be built into a build-type dir")
    SET(IFEM_LIBRARIES ${PROJECT_SOURCE_DIR}/../../${CMAKE_BUILD_TYPE}/lib/libIFEM.a)
  ENDIF(NOT IFEM_LIBRARIES)
ELSE(IFEM_INCLUDES)
  MESSAGE(STATUS "No in-tree libIFEM found, looking for system library")
  FIND_PATH(IFEM_INCLUDES
            NAMES
            SIMbase.h
            PATHS
            $ENV{HOME}/include
            PATH_SUFFIXES IFEM)
  FIND_LIBRARY(IFEM_LIBRARIES 
               IFEM
               PATHS
               $ENV{HOME}/lib)
ENDIF(IFEM_INCLUDES)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IFEM DEFAULT_MSG
                                  IFEM_INCLUDES IFEM_LIBRARIES)

MARK_AS_ADVANCED(IFEM_INCLUDES IFEM_LIBRARIES)
