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
elseif(EXISTS ${PROJECT_SOURCE_DIR}/Poisson)
  add_subdirectory(Poisson)
endif()

set(HAVE_ADVDIFF FALSE)
if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-AdvectionDiffusion)
  add_subdirectory(IFEM-AdvectionDiffusion)
  set(HAVE_ADVDIFF TRUE)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/AdvectionDiffusion)
  add_subdirectory(AdvectionDiffusion)
  set(HAVE_ADVDIFF TRUE)
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Darcy)
  add_subdirectory(IFEM-Darcy)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/Darcy)
  add_subdirectory(Darcy)
endif()

set(HAVE_BEAMEX FALSE)
set(HAVE_FINITEDEF FALSE)
if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity)
  set(ELASTICITY_SRC IFEM-Elasticity)
else()
  set(ELASTICITY_SRC Elasticity)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC})

  add_subdirectory(${ELASTICITY_SRC})
  add_subdirectory(${ELASTICITY_SRC}/Beam)
  add_subdirectory(${ELASTICITY_SRC}/Linear)

  if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/IFEM-BeamEx)
    set(BEAMEX_SRC ${ELASTICITY_SRC}/IFEM-BeamEx)
  else()
    set(BEAMEX_SRC ${ELASTICITY_SRC}/BeamEx)
  endif()
  if (EXISTS ${PROJECT_SOURCE_DIR}/${BEAMEX_SRC})
    add_subdirectory(${BEAMEX_SRC})
    add_subdirectory(${BEAMEX_SRC}/BeamSim)
    set(HAVE_BEAMEX TRUE)
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/IFEM-FiniteDeformation)
    set(FINITEDEF_SRC ${ELASTICITY_SRC}/IFEM-FiniteDeformation)
  else()
    set(FINITEDEF_SRC ${ELASTICITY_SRC}/FiniteDeformation)
  endif()
  if (EXISTS ${PROJECT_SOURCE_DIR}/${FINITEDEF_SRC})
    add_subdirectory(${FINITEDEF_SRC})
    add_subdirectory(${FINITEDEF_SRC}/Nonlinear)
    set(HAVE_FINITEDEF TRUE)
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-ThermoElasticity)
    add_subdirectory(IFEM-ThermoElasticity)
  elseif(EXISTS ${PROJECT_SOURCE_DIR}/ThermoElasticity)
    add_subdirectory(ThermoElasticity)
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-PoroElasticity)
    add_subdirectory(IFEM-PoroElasticity)
  elseif(EXISTS ${PROJECT_SOURCE_DIR}/PoroElasticity)
    add_subdirectory(PoroElasticity)
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-OpenFrac)
    add_subdirectory(IFEM-OpenFrac)
  elseif(EXISTS ${PROJECT_SOURCE_DIR}/FractureDynamics)
    add_subdirectory(FractureDynamics)
  endif()

endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Stokes)
  set(STOKES_SRC IFEM-Stokes)
else()
  set(STOKES_SRC Stokes)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/${STOKES_SRC})

  add_subdirectory(${PROJECT_SOURCE_DIR}/${STOKES_SRC})

  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-NavierStokes)
    set(NAVIERSTOKES_SRC IFEM-NavierStokes)
  else()
    set(NAVIERSTOKES_SRC NavierStokes)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/${NAVIERSTOKES_SRC})

    if (${HAVE_ADVDIFF})
      add_subdirectory(${NAVIERSTOKES_SRC})
    endif()
    if (${HAVE_BEAMEX} AND ${HAVE_FINITEDEF})

      if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-FSI)
        set(FSI_SRC IFEM-FSI)
      else()
        set(FSI_SRC FSI)
      endif()
      if(EXISTS ${PROJECT_SOURCE_DIR}/${FSI})
        add_subdirectory(${FSI})
      endif()

    endif()
  endif()
endif()

add_subdirectory(HDF5toVTx)

# Add 'check' target which builds all test applications, then executes the tests
add_check_target()
