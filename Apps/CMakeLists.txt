# This is a top-level CMakeLists.txt which will build all IFEM apps
# which exist as a local clone in sub-directories.
#
# NOTE: This file is _not_ used by the Jenkins build system on github.
# Therefore, make sure to test any changes locally before merging.

cmake_minimum_required(VERSION 3.22)

project(IFEM_Apps)

# Macro for including an app from two optional locations
macro(add_app NAME FOUND)
  set(${FOUND} FALSE)
  if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-${NAME})
    set(${FOUND} TRUE)
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/IFEM-${NAME}")
    add_subdirectory(IFEM-${NAME})
  elseif(EXISTS ${PROJECT_SOURCE_DIR}/${NAME})
    set(${FOUND} TRUE)
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/${NAME}")
    add_subdirectory(${NAME})
  endif(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-${NAME})
endmacro(add_app)

# Add local modules
set(IFEM_COMMON_APP_BUILD 1)

if(IFEM_AS_SUBMODULE)
  add_subdirectory(.. IFEM)
  set(IFEM_PATH ${PROJECT_SOURCE_DIR}/..)
else()
  get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
  set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../${BUILD_DIR})
  find_package(IFEM REQUIRED)
endif()
enable_testing()

add_app(Poisson HAVE_POISSON)
add_app(AdvectionDiffusion HAVE_ADVDIFF)
add_app(Darcy HAVE_DARCY)

set(HAVE_BEAMEX FALSE)
set(HAVE_FINITEDEF FALSE)
if(EXISTS ${PROJECT_SOURCE_DIR}/IFEM-Elasticity)
  set(ELASTICITY_SRC IFEM-Elasticity)
else()
  set(ELASTICITY_SRC Elasticity)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC})

  message(STATUS "Including ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}")
  add_subdirectory(${ELASTICITY_SRC})
  message(STATUS "Including ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/Beam")
  add_subdirectory(${ELASTICITY_SRC}/Beam)
  message(STATUS "Including ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/Shell")
  add_subdirectory(${ELASTICITY_SRC}/Shell)

  if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/IFEM-BeamGeometry)
    set(BEAMGEO_SRC ${ELASTICITY_SRC}/IFEM-BeamGeometry)
  else()
    set(BEAMGEO_SRC ${ELASTICITY_SRC}/BeamGeometry)
  endif()
  if (EXISTS ${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC})
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC}")
    add_subdirectory(${BEAMGEO_SRC})
  endif()

  message(STATUS "Including ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/Linear")
  add_subdirectory(${ELASTICITY_SRC}/Linear)

  if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/IFEM-BeamEx)
    set(BEAMEX_SRC ${ELASTICITY_SRC}/IFEM-BeamEx)
  else()
    set(BEAMEX_SRC ${ELASTICITY_SRC}/BeamEx)
  endif()
  if (EXISTS ${PROJECT_SOURCE_DIR}/${BEAMEX_SRC})
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/${BEAMEX_SRC}")
    add_subdirectory(${BEAMEX_SRC})
    add_subdirectory(${BEAMEX_SRC}/BeamSim)
    set(HAVE_BEAMEX TRUE)
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/IFEM-ShellEx)
    set(SHELLEX_SRC ${ELASTICITY_SRC}/IFEM-ShellEx)
  else()
    set(SHELLEX_SRC ${ELASTICITY_SRC}/ShellEx)
  endif()
  if (EXISTS ${PROJECT_SOURCE_DIR}/${SHELLEX_SRC})
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/${SHELLEX_SRC}")
    add_subdirectory(${SHELLEX_SRC})
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/${ELASTICITY_SRC}/IFEM-FiniteDeformation)
    set(FINITEDEF_SRC ${ELASTICITY_SRC}/IFEM-FiniteDeformation)
  else()
    set(FINITEDEF_SRC ${ELASTICITY_SRC}/FiniteDeformation)
  endif()
  if (EXISTS ${PROJECT_SOURCE_DIR}/${FINITEDEF_SRC})
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/${FINITEDEF_SRC}")
    add_subdirectory(${FINITEDEF_SRC})
    add_subdirectory(${FINITEDEF_SRC}/Nonlinear)
    set(HAVE_FINITEDEF TRUE)
  endif()

  add_app(ThermoElasticity HAVE_THERMOEL)
  add_app(PoroElasticity HAVE_POROEL)
  if(HAVE_POROEL)
    add_app(THM HAVE_THM)
  endif()
  add_app(OpenFrac HAVE_OPENFRAC)
  if(NOT HAVE_OPENFRAC AND EXISTS ${PROJECT_SOURCE_DIR}/FractureDynamics)
    message(STATUS "Including ${PROJECT_SOURCE_DIR}/FractureDynamics")
    add_subdirectory(FractureDynamics)
  endif()

endif()

add_app(Stokes HAVE_STOKES)
if(HAVE_STOKES)
  add_app(NavierStokes HAVE_NAVIERSTOKES)
  if(HAVE_NAVIERSTOKES)
    if (HAVE_ADVDIFF)
      add_app(Boussinesq HAVE_BOUSSINESQ)
    endif()
    if (HAVE_BEAMEX AND HAVE_FINITEDEF)
      add_app(FSI HAVE_FSI)
    endif()
  endif()
endif()

add_app(SIMRA-PostProc HAVE_POSTPROC)
add_app(HM HAVE_HM)
add_app(CoSTA HAVE_COSTA)
add_app(Burgers HAVE_BURGERS)

# Add 'check' target which builds all test applications,
# then executes all the tests
ifem_add_check_target()
