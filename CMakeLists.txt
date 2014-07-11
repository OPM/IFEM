PROJECT(IFEM)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(IFEM_VERSION_MAJOR 0)
SET(IFEM_VERSION_MINOR 9)
SET(IFEM_VERSION_PATCH 0)
SET(IFEM_VERSION "${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR}.${IFEM_VERSION_PATCH}")
SET(IFEM_ABI_VERSION ${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR})

# Generate header with version info
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/src/IFEM.h.in
               IFEM.h @ONLY)

# Generate doxy in build tree
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Doxyfile.in
               Doxyfile @ONLY)

# Required defines (well, PROFILE_LEVEL is not actually required, but...)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DReal=double -DepsZ=1.0e-12 -DPROFILE_LEVEL=3")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DINDEX_CHECK=2")
IF(VERBOSE_DEBUG GREATER 0)
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSP_DEBUG=${VERBOSE_DEBUG}")
ENDIF(VERBOSE_DEBUG GREATER 0)

ENABLE_TESTING()

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_SOURCE_DIR}/cmake/Modules)

INCLUDE(Apps/Common/cmake/UseMultiArch.cmake)

FIND_PACKAGE(IFEMDeps)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_BUILD_CXX_FLAGS}")

IF(NOT IFEM_USE_SYSTEM_TINYXML)
  SET(IFEM_INCLUDES ${IFEM_INCLUDES}
                    ${PROJECT_SOURCE_DIR}/3rdparty/tinyxml)
  FILE(GLOB_RECURSE TINYXML_SRCS ${PROJECT_SOURCE_DIR}/3rdparty/tinyxml/*.C)
ENDIF(NOT IFEM_USE_SYSTEM_TINYXML)

IF(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)
  SET(IFEM_INCLUDES ${IFEM_INCLUDES}
                    ${PROJECT_SOURCE_DIR}/src/ASM/LR)
ENDIF(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)


# Required include directories
SET(IFEM_INCLUDES
    ${IFEM_INCLUDES}
    ${IFEM_DEPINCLUDES}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/ASM
    ${PROJECT_SOURCE_DIR}/src/Eig
    ${PROJECT_SOURCE_DIR}/src/LinAlg
    ${PROJECT_SOURCE_DIR}/src/SIM
    ${PROJECT_SOURCE_DIR}/src/Utility
    ${PROJECT_SOURCE_DIR}/3rdparty/expreval
    ${PROJECT_BINARY_DIR}
)

INCLUDE_DIRECTORIES(${IFEM_INCLUDES})

SET(EXECUTABLE_OUTPUT_PATH bin)
SET(LIBRARY_OUTPUT_PATH lib)

IF(NOT CMAKE_INSTALL_DOCDIR)
  SET(CMAKE_INSTALL_DOCDIR share/doc/libifem1)
ENDIF(NOT CMAKE_INSTALL_DOCDIR)

IF(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-parentheses")
ENDIF(NOT WIN32)

# Generate build flags in build tree
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/cmake/Modules/IFEMFlags.cmake.in
               IFEMFlags.cmake @ONLY)

# Make the IFEM library
FILE(GLOB_RECURSE IFEM_SRCS ${PROJECT_SOURCE_DIR}/src/*.C
                            ${PROJECT_SOURCE_DIR}/src/*.f
                            ${PROJECT_SOURCE_DIR}/3rdparty/expreval/*.cpp)
string(REGEX REPLACE "${PROJECT_SOURCE_DIR}/src/ASM/ASMs.DInterpolate.C" "" IFEM_SRCS "${IFEM_SRCS}")
IF(NOT LRSPLINE_FOUND AND NOT LRSpline_FOUND)
  string(REGEX REPLACE "${PROJECT_SOURCE_DIR}/src/ASM/LR/[^;]*" "" IFEM_SRCS "${IFEM_SRCS}")
ENDIF()
ADD_LIBRARY(IFEM ${IFEM_SRCS} ${TINYXML_SRCS})
TARGET_LINK_LIBRARIES(IFEM ${IFEM_DEPLIBS})
SET_TARGET_PROPERTIES(IFEM PROPERTIES VERSION ${IFEM_VERSION}
                           SOVERSION ${IFEM_ABI_VERSION})

# For generating the doxy
ADD_CUSTOM_TARGET(lib_doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  COMMENT "Generating API documentation" VERBATIM)
if(NOT TARGET doc)
  add_custom_target(doc)
endif()
add_dependencies(doc lib_doc)

# 'install' target
# Have install force doc build
INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)

IF(WIN32)
  # TODO
ELSE(WIN32)
  # lib
  INSTALL(TARGETS IFEM DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT bin)

  # headers
  FILE(GLOB IFEM_HEADERS src/ASM/AlgEqSystem.h src/ASM/ASMbase.h
                         src/ASM/ASMenums.h src/ASM/ASMmxBase.h
                         src/ASM/ASMstruct.h src/ASM/ElmMats.h src/ASM/ElmNorm.h
                         src/ASM/Field.h src/ASM/Fields.h
                         src/ASM/FiniteElement.h src/ASM/GlbNorm.h
                         src/ASM/GlobalIntegral.h src/ASM/IntegrandBase.h
                         src/ASM/Integrand.h src/ASM/Lagrange.h
                         src/ASM/LocalIntegral.h src/ASM/SAMpatch.h
                         src/ASM/TimeDomain.h src/ASM/ASMs?D.h src/ASM/ASM?D.h
                         src/Integrands/*.h src/LinAlg/*.h src/SIM/*.h
                         src/Utility/*.h
                         ${CMAKE_BINARY_DIR}/IFEM.h)
  INSTALL(FILES ${IFEM_HEADERS}
    DESTINATION include/IFEM COMPONENT ${IFEM_DEV_COMPONENT})

  # cmake modules
  INSTALL(FILES cmake/Modules/FindIFEM.cmake DESTINATION lib/IFEM RENAME IFEMConfig.cmake)

  # documentation and license
  IF(INSTALL_DOXY)
    INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
    INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}
            COMPONENT doc
            PATTERN *.md5 EXCLUDE
            PATTERN *.map EXCLUDE)
    INSTALL(FILES doc/petsc.pdf doc/sim-input.pdf doc/Tutorial/GettingStarted.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})
  ENDIF(INSTALL_DOXY)
ENDIF(WIN32)
