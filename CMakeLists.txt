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
endif()

# Required include directories
set(IFEM_INCLUDES
    ${IFEM_INCLUDES}
    ${IFEM_DEPINCLUDES}
    ${PROJECT_SOURCE_DIR}/src/ASM
    ${PROJECT_SOURCE_DIR}/src/Eig
    ${PROJECT_SOURCE_DIR}/src/LinAlg
    ${PROJECT_SOURCE_DIR}/src/SIM
    ${PROJECT_SOURCE_DIR}/src/Utility
    ${PROJECT_SOURCE_DIR}/3rdparty
    ${PROJECT_SOURCE_DIR}/3rdparty/expreval
    ${PROJECT_SOURCE_DIR}/3rdparty/gtest/include
    ${PROJECT_BINARY_DIR}
)

include_directories(${IFEM_INCLUDES})

set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH lib)

# Generate build flags in build tree
configure_file(${PROJECT_SOURCE_DIR}/cmake/Modules/IFEMFlags.cmake.in
               IFEMFlags.cmake @ONLY)

# Make the IFEM library
file(GLOB IFEM_SRCS ${PROJECT_SOURCE_DIR}/src/IFEM.C
                    ${PROJECT_SOURCE_DIR}/src/ASM/*.C
                    ${PROJECT_SOURCE_DIR}/src/Eig/*.C
                    ${PROJECT_SOURCE_DIR}/src/LinAlg/*.C
                    ${PROJECT_SOURCE_DIR}/src/SIM/*.C
                    ${PROJECT_SOURCE_DIR}/src/Utility/*.C
                    ${PROJECT_SOURCE_DIR}/src/ASM/*.f
                    ${PROJECT_SOURCE_DIR}/src/Eig/*.f
                    ${PROJECT_SOURCE_DIR}/3rdparty/expreval/*.cpp)
string(REGEX REPLACE "${PROJECT_SOURCE_DIR}/src/ASM/ASMs.DInterpolate.C" "" IFEM_SRCS "${IFEM_SRCS}")

if(NOT PETSC_FOUND)
  list(REMOVE_ITEM IFEM_SRCS
                   ${PROJECT_SOURCE_DIR}/src/ASM/SAMpatchPETSc.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/PETScMatrix.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/PETScPCPerm.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/PETScPCProd.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/PETScPCScale.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/PETScSchurPC.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/PETScSolParams.C)
endif()

if(NOT ISTL_FOUND)
  list(REMOVE_ITEM IFEM_SRCS
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/ISTLMatrix.C
                   ${PROJECT_SOURCE_DIR}/src/LinAlg/ISTLSolParams.C)
endif()

if(LRSPLINE_FOUND OR LRSpline_FOUND)
  file(GLOB LR_SRCS ${PROJECT_SOURCE_DIR}/src/ASM/LR/*.C)
  list(APPEND IFEM_SRCS ${LR_SRCS})
endif()
list(APPEND CHECK_SOURCES ${IFEM_SRCS})
add_library(IFEM ${IFEM_SRCS} ${TINYXML_SRCS})
target_link_libraries(IFEM ${IFEM_DEPLIBS})
set_target_properties(IFEM PROPERTIES VERSION ${IFEM_VERSION}
                           SOVERSION ${IFEM_ABI_VERSION})

set(IFEM_PATH ${PROJECT_SOURCE_DIR})
if(NOT IFEM_SUBTREE_BUILD AND NOT IFEM_AS_SUBMODULE)
  set(EXTRA_DOXY_PATHS ${PROJECT_SOURCE_DIR}/doc)
endif()

# Generate doxy in build tree
configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
               Doxyfile @ONLY)

# For generating the doxy
add_custom_target(lib_doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  COMMENT "Generating API documentation" VERBATIM)
if(NOT TARGET doc AND NOT IFEM_AS_SUBMODULE)
  add_custom_target(doc)
  add_dependencies(doc lib_doc)
endif()

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
                         src/ASM/ASM*truct.h src/ASM/*Mats.h src/ASM/ElmNorm.h
                         src/ASM/Field.h src/ASM/Fields.h src/ASM/GlbForceVec.h
                         src/ASM/FiniteElement.h src/ASM/GlbNorm.h
                         src/ASM/GlobalIntegral.h src/ASM/IntegrandBase.h
                         src/ASM/ImmersedBoundaries.h src/ASM/Interface.h
                         src/ASM/Integrand.h src/ASM/Lagrange.h
                         src/ASM/LocalIntegral.h src/ASM/SAMpatch.h
                         src/ASM/TimeDomain.h src/ASM/ASMs?D.h src/ASM/ASM?D.h
                         src/ASM/ASMs?DLag.h
                         src/ASM/DomainDecomposition.h src/ASM/ItgPoint.h
                         src/ASM/ReactionsOnly.h
                         src/LinAlg/*.h src/SIM/*.h src/Utility/*.h
                         3rdparty/*.h
                         ${CMAKE_BINARY_DIR}/IFEM.h)
  if(NOT IFEM_USE_SYSTEM_TINYXML)
    list(APPEND IFEM_HEADERS 3rdparty/tinyxml/tinyxml.h)
  endif()
  install(FILES ${IFEM_HEADERS}
          DESTINATION include/IFEM COMPONENT ${IFEM_DEV_COMPONENT})

  set(IFEM_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include/IFEM ${IFEM_DEPINCLUDES})
  if(LRSpline_FOUND)
    file(GLOB IFEM_LR_HEADERS src/ASM/LR/*.h)
    install(FILES ${IFEM_LR_HEADERS} DESTINATION include/IFEM/LR)
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
                cmake/Scripts/IFEMDoxy.cmake
                cmake/Scripts/IFEMTesting.cmake DESTINATION lib/IFEM)

  # source code
  install(FILES src/IFEM-test.C DESTINATION lib/IFEM/src)

  # documentation and license
  if(IFEM_INSTALL_DOXY)
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
  set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
endif()
