cmake_minimum_required(VERSION 3.22)

project(IFEMAppCommon LANGUAGES C CXX)

# Todo: Remove list when we can drop support for cmake 3.22
set(APPCOMMON_PUBLIC_HEADERS
  CompatibleOperators.h
  EqualOrderOperators.h
  ISolver.h
  MeshUtils.h
  MultiPatchModelGenerator.h
  PiolaOperators.h
  SIMconfigure.h
  SIMCoupled.h
  SIMCoupledSI.h
  SIMExplicitLMM.h
  SIMExplicitRKE.h
  SIMExplicitRK.h
  SIMImplicitLMM.h
  SIMMultiPatchModelGen.h
  SIMNodalConstraint.h
  SIMSemi3D.h
  SIMSolverAdap.h
  SIMSolver.h
  SIMSolverTS.h
  StabilizationUtils.h
  TimeIntUtils.h
)

ifem_add_library(
  NAME
    IFEMAppCommon
  SOURCES
    CompatibleOperators.C
    EqualOrderOperators.C
    MeshUtils.C
    MultiPatchModelGenerator.C
    PiolaOperators.C
    SIMMultiPatchModelGen.C
    SIMNodalConstraint.C
    SIMSemi3D.C
    StabilizationUtils.C
    TimeIntUtils.C
  HEADERS
    ${APPCOMMON_PUBLIC_HEADERS}
  LIBRARIES
    IFEM
)

# Unit tests
ifem_add_test_app(
  NAME
    AppCommon
  SOURCES
    Test/TestCompatibleOperators.C
    Test/TestEqualOrderOperators.C
    Test/TestMeshUtils.C
    Test/TestMultiPatchModelGenerator.C
    Test/TestSIMCoupled.C
    Test/TestSIMNodalConstraint.C
    Test/TestStabilizationUtils.C
    Test/TestTimeIntUtils.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test
  LIBRARIES
    IFEMAppCommon
)

if(LRSpline_FOUND AND TARGET AppCommon-test)
  target_sources(
    AppCommon-test PRIVATE
    Test/LR/TestMultiPatchLRRefine.C
    Test/LR/TestSIMNodalConstraint.C
  )
endif()

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)

# Installation
include(GNUInstallDirs)

target_include_directories(IFEMAppCommon INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/IFEM>)

if(CMAKE_VERSION GREATER_EQUAL 3.23)
  install(TARGETS IFEMAppCommon
          EXPORT IFEM
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
          FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/IFEM)
else()
  install(TARGETS IFEMAppCommon
          EXPORT IFEM
          DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(
    FILES
      ${APPCOMMON_PUBLIC_HEADERS}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/IFEM
  )
endif()

install(FILES scripts/regtest.sh
              scripts/iotest.sh
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/IFEM/scripts
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE)
