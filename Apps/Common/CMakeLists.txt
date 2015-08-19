PROJECT(CommonIFEM)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Required packages
IF(NOT IFEM_CONFIGURED)
  FIND_PACKAGE(IFEM REQUIRED)
ENDIF(NOT IFEM_CONFIGURED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR})

IF(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(NOT WIN32)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

INCLUDE(cmake/UseMultiArch.cmake)

# Common application sources
ADD_LIBRARY(CommonIFEM MeshUtils.C
                       ResidualOperators.C
                       Spalding.C
                       SAWallLaw.C
                       StabilizationUtils.C
                       TimeIntUtils.C
                       WeakOperators.C)
TARGET_LINK_LIBRARIES(CommonIFEM ${IFEM_LIBRARIES})
SET_TARGET_PROPERTIES(CommonIFEM PROPERTIES
                      VERSION ${IFEM_VERSION}
                      SOVERSION ${IFEM_ABI_VERSION})
INSTALL(TARGETS CommonIFEM DESTINATION ${CMAKE_INSTALL_LIBDIR})
FILE(GLOB HEADERS *.h)
INSTALL(FILES ${HEADERS} DESTINATION include/IFEM)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  AppCommon
                  CommonIFEM ${IFEM_LIBRARIES})
set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
