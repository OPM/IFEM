project(IFEM_Apps)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(IFEM_PATH ${PROJECT_SOURCE_DIR}/..)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${IFEM_PATH}/cmake/Modules)
set(IFEM_COMMON_APP_BUILD 1)

find_package(IFEM REQUIRED)

enable_testing()

include_directories(${IFEM_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
if(NOT IFEM_AS_SUBMODULE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROFILE_LEVEL=2")
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Poisson)
  add_subdirectory(IFEM-Poisson)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-AdvectionDiffusion)
        add_subdirectory(IFEM-AdvectionDiffusion)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Darcy)
        add_subdirectory(IFEM-Darcy)
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity)
        add_subdirectory(IFEM-Elasticity)
        add_subdirectory(IFEM-Elasticity/Beam)
        add_subdirectory(IFEM-Elasticity/Linear)
        if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity/IFEM-BeamEx)
                add_subdirectory(IFEM-Elasticity/IFEM-BeamEx)
                add_subdirectory(IFEM-Elasticity/IFEM-BeamEx/BeamSim)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity/IFEM-FiniteDeformation)
          add_subdirectory(IFEM-Elasticity/IFEM-FiniteDeformation)
          add_subdirectory(IFEM-Elasticity/IFEM-FiniteDeformation/Nonlinear)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-ThermoElasticity)
          add_subdirectory(IFEM-ThermoElasticity)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-PoroElasticity)
    add_subdirectory(IFEM-PoroElasticity)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-OpenFrac)
    add_subdirectory(IFEM-OpenFrac)
  endif()
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Stokes)
  add_subdirectory(${PROJECT_SOURCE_DIR}/IFEM-Stokes)
  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-NavierStokes)
    if (EXISTS ${PROJECT_SOURCE_DIR}/IFEM-AdvectionDiffusion)
      add_subdirectory(IFEM-NavierStokes)
    endif()
    if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-FSI AND
        EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity/IFEM-BeamEx AND
        EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity/IFEM-FiniteDeformation)
      add_subdirectory(IFEM-FSI)
    endif()
  endif()
endif()

add_subdirectory(HDF5toVTx)

# Add 'check' target which builds all test applications, then executes the tests
add_check_target()
