PROJECT(IFEM_Apps)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../cmake/Modules)

FIND_PACKAGE(IFEM REQUIRED)
ENABLE_TESTING()

ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(AdvectionDiffusion)
ADD_SUBDIRECTORY(FSWallDistance)
ADD_SUBDIRECTORY(Elasticity/Linear)
ADD_SUBDIRECTORY(Elasticity/FiniteDeformation)
ADD_SUBDIRECTORY(HDF5toVTx)
ADD_SUBDIRECTORY(Poisson)
ADD_SUBDIRECTORY(Stokes)

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Common/scripts/regtest.sh.in
               regtest.sh)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Common/scripts/doregtests.sh.in
  doregtests.sh @ONLY)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Common/scripts/check-patchseries.sh.in
               check-patchseries.sh)
