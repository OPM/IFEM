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

if(EXISTS ${PROJECT_SOURCE_DIR}/Poisson)
  add_subdirectory(Poisson)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/AdvectionDiffusion)
  add_subdirectory(AdvectionDiffusion)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/Darcy)
  add_subdirectory(Darcy)
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/Elasticity)
  add_subdirectory(Elasticity)
  add_subdirectory(Elasticity/Beam)
  add_subdirectory(Elasticity/Linear)
  if(EXISTS ${PROJECT_SOURCE_DIR}/Elasticity/BeamEx)
    add_subdirectory(Elasticity/BeamEx)
    add_subdirectory(Elasticity/BeamSim)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/Elasticity/FiniteDeformation)
    add_subdirectory(Elasticity/FiniteDeformation)
    add_subdirectory(Elasticity/Nonlinear)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/ThermoElasticity)
    add_subdirectory(ThermoElasticity)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/PoroElasticity)
    add_subdirectory(PoroElasticity)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/FractureDynamics)
    add_subdirectory(FractureDynamics)
  endif()
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/Stokes)
  add_subdirectory(${PROJECT_SOURCE_DIR}/Stokes)
  if(EXISTS ${PROJECT_SOURCE_DIR}/NavierStokes)
    if (EXISTS ${PROJECT_SOURCE_DIR}/AdvectionDiffusion)
      add_subdirectory(NavierStokes)
    endif()
    if(EXISTS ${PROJECT_SOURCE_DIR}/FSI AND
       EXISTS ${PROJECT_SOURCE_DIR}/Elasticity/BeamEx AND
       EXISTS ${PROJECT_SOURCE_DIR}/Elasticity/FiniteDeformation)
      add_subdirectory(FSI)
    endif()
  endif()
endif()

add_subdirectory(HDF5toVTx)

# Add 'check' target which builds all test applications, then executes the tests
add_check_target()
