PROJECT(IFEM_Apps)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../cmake/Modules)

FIND_PACKAGE(IFEM REQUIRED)
ENABLE_TESTING()

include_directories(${IFEM_INCLUDES})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
SET(EXECUTABLE_OUTPUT_PATH bin)

set(IFEM_COMMON_APP_BUILD 1)
set(IFEM_PATH ${PROJECT_SOURCE_DIR}/..)
ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(Poisson)
ADD_SUBDIRECTORY(AdvectionDiffusion)
ADD_SUBDIRECTORY(Elasticity/Linear)
ADD_SUBDIRECTORY(Elasticity/Dynamic)
ADD_SUBDIRECTORY(Elasticity/FiniteDeformation)
ADD_SUBDIRECTORY(FSWallDistance)
ADD_SUBDIRECTORY(Stokes)
ADD_SUBDIRECTORY(ThermoElasticity)
ADD_SUBDIRECTORY(HDF5toVTx)

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Common/scripts/regtest.sh.in
               regtest.sh)

add_check_target()
