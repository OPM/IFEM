project(CommonIFEM)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Required packages
if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR})

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# Common application sources
add_library(CommonIFEM STATIC MeshUtils.C
                       ResidualOperators.C
                       Spalding.C
                       SAWallLaw.C
                       SIMNodalConstraint.C
                       SIMSemi3D.C
                       StabilizationUtils.C
                       TimeIntUtils.C
                       WeakOperators.C)
target_link_libraries(CommonIFEM ${IFEM_LIBRARIES})

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  AppCommon
                  CommonIFEM ${IFEM_LIBRARIES})
set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
