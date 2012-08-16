IF(NOT DEFINED FORCE_SYSTEM_IFEM OR NOT "${FORCE_SYSTEM_IFEM}")
  FIND_PATH(IFEM_PATH
    NAMES IFEM.h
    PATHS ${PROJECT_SOURCE_DIR}/../../src/
          ${PROJECT_SOURCE_DIR}/../../../src/)
ENDIF(NOT DEFINED FORCE_SYSTEM_IFEM OR NOT "${FORCE_SYSTEM_IFEM}")

FIND_PACKAGE(IFEMDeps)

IF(IFEM_PATH)
  # Build is in-tree
  MESSAGE(STATUS "Using in-tree libIFEM")

  SET(IFEM_INCLUDES ${IFEM_PATH}
                    ${IFEM_PATH}/ASM
                    ${IFEM_PATH}/Eig
                    ${IFEM_PATH}/LinAlg
                    ${IFEM_PATH}/SIM
                    ${IFEM_PATH}/Utility)

  FIND_LIBRARY(IFEM_LIBRARIES
    NAMES IFEM
    PATHS ${IFEM_PATH}/../${IFEM_BUILD_TYPE}/lib
          ${IFEM_PATH}/../lib
    NO_DEFAULT_PATH)
  IF(NOT IFEM_LIBRARIES)
    MESSAGE(WARNING "Could not find the in-tree libIFEM library, "
      "we assume it will be built into a build-type dir")
  ENDIF(NOT IFEM_LIBRARIES)

  IF(NOT IFEM_USE_SYSTEM_TINYXML)
    SET(IFEM_INCLUDES ${IFEM_INCLUDES}
                      ${IFEM_PATH}/../3rdparty/tinyxml)
  ENDIF(NOT IFEM_USE_SYSTEM_TINYXML)
ELSE(IFEM_PATH)
  IF(NOT DEFINED FORCE_SYSTEM_IFEM OR NOT "${FORCE_SYSTEM_IFEM}")
    MESSAGE(STATUS "No in-tree libIFEM found, looking for system library")
  ENDIF(NOT DEFINED FORCE_SYSTEM_IFEM OR NOT "${FORCE_SYSTEM_IFEM}")

  FIND_PATH(IFEM_INCLUDES
    NAMES SIMbase.h
    PATHS $ENV{HOME}/include
    PATH_SUFFIXES IFEM)

  FIND_LIBRARY(IFEM_LIBRARIES
    NAMES IFEM
    PATHS $ENV{HOME}/lib)
ENDIF(IFEM_PATH)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IFEM DEFAULT_MSG
                                  IFEM_LIBRARIES)

SET(IFEM_LIBRARIES ${IFEM_LIBRARIES} ${IFEM_DEPLIBS})
SET(IFEM_INCLUDES ${IFEM_INCLUDES} ${IFEM_DEPINCLUDES})

