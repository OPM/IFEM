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

# Required include directories
SET(INCLUDES
  ${IFEM_DEPINCLUDES}
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/src/ASM
  ${PROJECT_SOURCE_DIR}/src/Eig
  ${PROJECT_SOURCE_DIR}/src/LinAlg
  ${PROJECT_SOURCE_DIR}/src/SIM
  ${PROJECT_SOURCE_DIR}/src/Utility
  ${PROJECT_SOURCE_DIR}/3rdparty/expreval
  ${CMAKE_BINARY_DIR}
)

IF(NOT IFEM_USE_SYSTEM_TINYXML)
  SET(INCLUDES ${INCLUDES}
               ${PROJECT_SOURCE_DIR}/3rdparty/tinyxml)
  FILE(GLOB_RECURSE TINYXML_SRCS ${PROJECT_SOURCE_DIR}/3rdparty/tinyxml/*.C)
ENDIF(NOT IFEM_USE_SYSTEM_TINYXML)

IF(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)
  SET(INCLUDES ${INCLUDES}
               ${PROJECT_SOURCE_DIR}/src/ASM/LR)
ENDIF(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)

INCLUDE_DIRECTORIES(${INCLUDES})

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

# Make the IFEM library
FILE(GLOB_RECURSE IFEM_SRCS ${PROJECT_SOURCE_DIR}/src/*.[Cf]
                            ${PROJECT_SOURCE_DIR}/3rdparty/expreval/*.cpp)
IF(NOT(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS))
  string(REGEX REPLACE "${PROJECT_SOURCE_DIR}/src/ASM/LR/[^;]*" "" IFEM_SRCS "${IFEM_SRCS}")
ENDIF(NOT(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS))
ADD_LIBRARY(IFEM ${IFEM_SRCS} ${TINYXML_SRCS})
TARGET_LINK_LIBRARIES(IFEM ${IFEM_DEPLIBS})
SET_TARGET_PROPERTIES(IFEM PROPERTIES VERSION ${IFEM_VERSION}
                           SOVERSION ${IFEM_ABI_VERSION})

# For generating the doxy
ADD_CUSTOM_TARGET(doc doxygen ${CMAKE_BINARY_DIR}/Doxyfile
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  COMMENT "Generating API documentation" VERBATIM)

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
                         src/ASM/TimeDomain.h
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
    INSTALL(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}
            COMPONENT doc
            PATTERN *.md5 EXCLUDE
            PATTERN *.map EXCLUDE)
  ENDIF(INSTALL_DOXY)
ENDIF(WIN32)
