project(IFEM)

cmake_minimum_required(VERSION 2.6)

set(IFEM_VERSION_MAJOR 0)
set(IFEM_VERSION_MINOR 9)
set(IFEM_VERSION_PATCH 0)
set(IFEM_VERSION "${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR}.${IFEM_VERSION_PATCH}")
set(IFEM_ABI_VERSION ${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR})

# Generate header with version info
configure_file(${PROJECT_SOURCE_DIR}/src/IFEM.h.in
               IFEM.h @ONLY)

# Generate doxy in build tree
configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
               Doxyfile @ONLY)

# Required defines (well, PROFILE_LEVEL is not actually required, but...)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DReal=double -DepsZ=1.0e-12 -DPROFILE_LEVEL=3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DINDEX_CHECK=2")
if(VERBOSE_DEBUG GREATER 0)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSP_DEBUG=${VERBOSE_DEBUG}")
endif()

set(IFEM_LIBRARY_BUILD 1)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_SOURCE_DIR}/cmake/Modules
                      ${CMAKE_SOURCE_DIR}/cmake/Scripts)

include(cmake/Scripts/UseMultiArch.cmake)

if(NOT IFEM_AS_SUBMODULE)
  if(NOT IFEM_PATH)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR})
  endif()
  find_package(IFEMDeps)
endif()

add_definitions(${IFEM_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_BUILD_CXX_FLAGS}")

if(NOT IFEM_USE_SYSTEM_TINYXML)
  set(IFEM_INCLUDES ${IFEM_INCLUDES}
                    ${PROJECT_SOURCE_DIR}/3rdparty/tinyxml)
  file(GLOB_RECURSE TINYXML_SRCS ${PROJECT_SOURCE_DIR}/3rdparty/tinyxml/*.C)
endif()

if(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)
  set(IFEM_INCLUDES ${IFEM_INCLUDES}
                    ${PROJECT_SOURCE_DIR}/src/ASM/LR)
ENDIF()

# Required include directories
set(IFEM_INCLUDES
    ${IFEM_INCLUDES}
    ${IFEM_DEPINCLUDES}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/ASM
    ${PROJECT_SOURCE_DIR}/src/Eig
    ${PROJECT_SOURCE_DIR}/src/LinAlg
    ${PROJECT_SOURCE_DIR}/src/SIM
    ${PROJECT_SOURCE_DIR}/src/Utility
    ${PROJECT_SOURCE_DIR}/3rdparty/expreval
    ${PROJECT_SOURCE_DIR}/3rdparty/gtest/include
    ${PROJECT_BINARY_DIR}
)

include_directories(${IFEM_INCLUDES})

set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH lib)

if(NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR share/doc/libifem1)
endif()

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-parentheses")
endif()

# Generate build flags in build tree
configure_file(${PROJECT_SOURCE_DIR}/cmake/Modules/IFEMFlags.cmake.in
               IFEMFlags.cmake @ONLY)

# Make the IFEM library
set(IFEM_ASM_SRCS src/ASM/AlgEqSystem.C
                  src/ASM/ASM1D.C
                  src/ASM/ASM2D.C
                  src/ASM/ASM3D.C
                  src/ASM/ASMbase.C
                  src/ASM/ASMmxBase.C
                  src/ASM/ASMs1D.C
                  src/ASM/ASMs1DC1.C
                  src/ASM/ASMs1DLag.C
                  src/ASM/ASMs1DSpec.C
                  src/ASM/ASMs2D.C
                  src/ASM/ASMs2DC1.C
                  src/ASM/ASMs2DIB.C
                  src/ASM/ASMs2DLag.C
                  src/ASM/ASMs2Dmx.C
                  src/ASM/ASMs2DmxLag.C
                  src/ASM/ASMs2Drecovery.C
                  src/ASM/ASMs2DSpec.C
                  src/ASM/ASMs3D.C
                  src/ASM/ASMs3DLag.C
                  src/ASM/ASMs3Dmx.C
                  src/ASM/ASMs3DmxLag.C
                  src/ASM/ASMs3Drecovery.C
                  src/ASM/ASMs3DSpec.C
                  src/ASM/ASMstruct.C
                  src/ASM/ASMunstruct.C
                  src/ASM/BDFMats.C
                  src/ASM/ElmMats.C
                  src/ASM/Field.C
                  src/ASM/Fields.C
                  src/ASM/GlbForceVec.C
                  src/ASM/GlbL2projector.C
                  src/ASM/GlbNorm.C
                  src/ASM/HHTMats.C
                  src/ASM/IBGeometries.C
                  src/ASM/ImmersedBoundaries.C
                  src/ASM/IntegrandBase.C
                  src/ASM/Lagrange.C
                  src/ASM/LagrangeField2D.C
                  src/ASM/LagrangeField3D.C
                  src/ASM/LagrangeFields2D.C
                  src/ASM/LagrangeFields3D.C
                  src/ASM/NewmarkMats.C
                  src/ASM/SAMpatch.C
                  src/ASM/SAMpatchPara.C
                  src/ASM/SplineField2D.C
                  src/ASM/SplineField3D.C
                  src/ASM/SplineFields2D.C
                  src/ASM/SplineFields3D.C
                  src/ASM/SplineInterpolator.C
                  src/ASM/dslbln.f
                  src/ASM/dsolv3.f)

set(IFEM_EIG_SRCS src/Eig/EigSolver.C
                   src/Eig/eig_drv1.f
                   src/Eig/eig_drv2.f
                   src/Eig/eig_drv3.f
                   src/Eig/eig_drv4.f
                   src/Eig/eig_drv5.f
                   src/Eig/eig_drv6.f)

set(IFEM_LINALG_SRCS src/LinAlg/DenseMatrix.C
                     src/LinAlg/LinAlgInit.C
                     src/LinAlg/LinSolParams.C
                     src/LinAlg/matrix.C
                     src/LinAlg/MatVec.C
                     src/LinAlg/PCPerm.C
                     src/LinAlg/PCProd.C
                     src/LinAlg/PCScale.C
                     src/LinAlg/PETScBlockMatrix.C
                     src/LinAlg/PETScMatrix.C
                     src/LinAlg/ProcessAdm.C
                     src/LinAlg/SAM.C
                     src/LinAlg/SparseMatrix.C
                     src/LinAlg/SPRMatrix.C
                     src/LinAlg/SystemMatrix.C)

set(IFEM_SIM_SRCS src/SIM/AdaptiveSIM.C
                   src/SIM/EigenModeSIM.C
                   src/SIM/ForceIntegrator.C
                   src/SIM/GenAlphaSIM.C
                   src/SIM/HHTSIM.C
                   src/SIM/InitialConditionHandler.C
                   src/SIM/MultiStepSIM.C
                   src/SIM/NewmarkNLSIM.C
                   src/SIM/NewmarkSIM.C
                   src/SIM/NonLinSIM.C
                   src/SIM/SIM1D.C
                   src/SIM/SIM2D.C
                   src/SIM/SIM3D.C
                   src/SIM/SIMbase.C
                   src/SIM/SIMdependency.C
                   src/SIM/SIMgeneric.C
                   src/SIM/SIMinput.C
                   src/SIM/SIMoptions.C
                   src/SIM/SIMoutput.C
                   src/SIM/TimeStep.C)

set(IFEM_UTIL_SRCS src/Utility/AnaSol.C
                   src/Utility/BDF.C
                   src/Utility/ControlFIFO.C
                   src/Utility/CoordinateMapping.C
                   src/Utility/DataExporter.C
                   src/Utility/ElementBlock.C
                   src/Utility/ExprFunctions.C
                   src/Utility/FieldFunctions.C
                   src/Utility/Function.C
                   src/Utility/Functions.C
                   src/Utility/GaussQuadrature.C
                   src/Utility/HDF5Writer.C
                   src/Utility/LagrangeInterpolator.C
                   src/Utility/Legendre.C
                   src/Utility/LogStream.C
                   src/Utility/Math.C
                   src/Utility/MPC.C
                   src/Utility/NodeVecFunc.C
                   src/Utility/Profiler.C
                   src/Utility/ScopedLogger.C
                   src/Utility/SplineUtils.C
                   src/Utility/StringUtils.C
                   src/Utility/Tensor4.C
                   src/Utility/Tensor.C
                   src/Utility/ThreadGroups.C
                   src/Utility/Utilities.C
                   src/Utility/Vec3Oper.C
                   src/Utility/VTF.C
                   src/Utility/XMLWriter.C)

set(THIRDPARTY_SRCS 3rdparty/expreval/except.cpp
                    3rdparty/expreval/expr.cpp
                    3rdparty/expreval/func.cpp
                    3rdparty/expreval/funclist.cpp
                    3rdparty/expreval/node.cpp
                    3rdparty/expreval/parser.cpp
                    3rdparty/expreval/vallist.cpp)

if(LRSPLINE_FOUND OR LRSpline_FOUND)
  list(APPEND IFEM_ASM_SRCS src/ASM/LR/ASMLRSpline.C
                            src/ASM/LR/ASMu2D.C
                            src/ASM/LR/ASMu2DIB.C
                            src/ASM/LR/ASMu2Dmx.C
                            src/ASM/LR/ASMu2Dmxrecovery.C
                            src/ASM/LR/ASMu2Drecovery.C
                            src/ASM/LR/ASMu3D.C)
endif()
add_library(IFEM src/IFEM.C ${IFEM_ASM_SRCS} ${IFEM_EIG_SRCS}
                            ${IFEM_LINALG_SRCS} ${IFEM_UTIL_SRCS}
                            ${IFEM_SIM_SRCS} ${THIRDPARTY_SRCS} ${TINYXML_SRCS})
target_link_libraries(IFEM ${IFEM_DEPLIBS})
set_target_properties(IFEM PROPERTIES VERSION ${IFEM_VERSION}
                           SOVERSION ${IFEM_ABI_VERSION})

# For generating the doxy
add_custom_target(lib_doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  COMMENT "Generating API documentation" VERBATIM)
if(NOT TARGET doc)
  add_custom_target(doc)
endif()
add_dependencies(doc lib_doc)

if(NOT IFEM_AS_SUBMODULE AND NOT IFEM_COMMON_APP_BUILD)
  set(IFEM_INTREE_BUILD ON)
  set(IFEM_LIBRARIES IFEM ${IFEM_DEPLIBS})
  set(IFEM_PATH ${PROJECT_SOURCE_DIR})
endif()

enable_testing()
include(cmake/Scripts/IFEMTesting.cmake)

add_subdirectory(Apps/Common)

# Unit tests
if(IFEM_AS_SUBMODULE OR IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  add_check_target()
endif()

if(WIN32)
  # TODO
else()
  # lib
  install(TARGETS IFEM DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT bin)

  # headers
  file(GLOB IFEM_HEADERS src/ASM/AlgEqSystem.h src/ASM/ASMbase.h
                         src/ASM/ASMenums.h src/ASM/ASMmxBase.h
                         src/ASM/ASMs?Dmx.h
                         src/ASM/ASMstruct.h src/ASM/*Mats.h src/ASM/ElmNorm.h
                         src/ASM/Field.h src/ASM/Fields.h src/ASM/GlbForceVec.h
                         src/ASM/FiniteElement.h src/ASM/GlbNorm.h
                         src/ASM/GlobalIntegral.h src/ASM/IntegrandBase.h
                         src/ASM/ImmersedBoundaries.h
                         src/ASM/Integrand.h src/ASM/Lagrange.h
                         src/ASM/LocalIntegral.h src/ASM/SAMpatch.h
                         src/ASM/TimeDomain.h src/ASM/ASMs?D.h src/ASM/ASM?D.h
                         src/LinAlg/*.h src/SIM/*.h
                         src/Utility/*.h
                         ${CMAKE_BINARY_DIR}/IFEM.h)
  install(FILES ${IFEM_HEADERS}
          DESTINATION include/IFEM COMPONENT ${IFEM_DEV_COMPONENT})

  set(IFEM_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include/IFEM ${IFEM_DEPINCLUDES})
  if(LRSpline_FOUND)
    file(GLOB IFEM_LR_HEADERS src/ASM/LR/*.h)
    install(FILES ${IFEM_LR_HEADERS} DESTINATION include/IFEM/LR)
    install(FILES src/ASM/ASMunstruct.h DESTINATION include/IFEM)
    list(APPEND IFEM_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include/IFEM/LR)
  endif()

  # cmake modules
  string(REPLACE ";" " " IFEM_CXX_FLAGS ${IFEM_CXX_FLAGS})
  string(REPLACE ";" " " IFEM_DEPLIBS "${IFEM_DEPLIBS}")
  string(REPLACE ";" " " IFEM_DEFINITIONS "${IFEM_DEFINITIONS}")
  configure_file(cmake/Templates/IFEMConfigVersion.cmake.in IFEMConfigVersion.cmake @ONLY)
  configure_file(cmake/Templates/IFEMConfig.cmake.in IFEMConfig.cmake @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/IFEMConfigVersion.cmake
                ${CMAKE_BINARY_DIR}/IFEMConfig.cmake
                cmake/Scripts/CheckCommits.cmake
                cmake/Scripts/DownloadGTest.cmake
                cmake/Scripts/UseMultiArch.cmake
                cmake/Scripts/IFEMTesting.cmake DESTINATION lib/IFEM)

  # source code
  install(FILES src/IFEM-test.C DESTINATION lib/IFEM/src)

  # documentation and license
  option(INSTALL_DOXY "Install doxy with installation target?" ON)
  if(INSTALL_DOXY)
    install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
    install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}
            COMPONENT doc
            PATTERN *.md5 EXCLUDE
            PATTERN *.map EXCLUDE)
    install(FILES doc/petsc.pdf doc/sim-input.pdf doc/Tutorials/GettingStarted.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})
  endif()
endif()

if(IFEM_COMMON_APP_BUILD OR IFEM_AS_SUBMODULE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
endif()
