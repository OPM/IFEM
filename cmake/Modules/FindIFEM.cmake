# IFEM uses C++ and Fortran
enable_language(CXX)
enable_language(Fortran)

# Custom profiles

if(NOT IFEM_BUILD_TYPE)
  set(IFEM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(${CMAKE_BUILD_TYPE} MATCHES "Nopt")
  set(CMAKE_BUILD_TYPE Debug)
elseif(${CMAKE_BUILD_TYPE} MATCHES "Debug-MPI")
  set(CMAKE_BUILD_TYPE Debug)
elseif(${CMAKE_BUILD_TYPE} MATCHES "Nomp")
  set(CMAKE_BUILD_TYPE Release)
elseif(${CMAKE_BUILD_TYPE} MATCHES "Release-MPI")
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT IFEM_AS_SUBMODULE)
  find_path(IFEM_PATH
    NAMES ${IFEM_BUILD_TYPE}/IFEM.h
    PATHS ${PROJECT_SOURCE_DIR}/..
          ${PROJECT_SOURCE_DIR}/../..
          ${PROJECT_SOURCE_DIR}/../../..
          ${PROJECT_SOURCE_DIR}/../../../..
          NO_DEFAULT_PATHS)
    set(IFEM_H_PATH ${IFEM_PATH}/${IFEM_BUILD_TYPE})
  if(NOT IFEM_PATH)
    find_path(IFEM_PATH
      NAMES IFEM.h
      PATHS ${PROJECT_SOURCE_DIR}/..
            ${PROJECT_SOURCE_DIR}/../..
            ${PROJECT_SOURCE_DIR}/../../..
            ${PROJECT_SOURCE_DIR}/../../../..
            NO_DEFAULT_PATHS)
    set(IFEM_H_PATH ${IFEM_PATH})
  endif(NOT IFEM_PATH)
endif()

if(IFEM_AS_SUBMODULE)
  # Build wants IFEM as a submodule
  message(STATUS "Building IFEM as a submodule")
  if(CMAKE_CROSSCOMPILING)
    set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE)
  endif()
  find_path(IFEM_PATH NAMES src/IFEM.h.in
                      PATHS ${PROJECT_SOURCE_DIR}/..
                            ${PROJECT_SOURCE_DIR}/../..
                            ${PROJECT_SOURCE_DIR}/../../..
                            ${PROJECT_SOURCE_DIR}/../../../.. NO_DEFAULT_PATHS)
  if(NOT IFEM_PATH)
    message(FATAL_ERROR "IFEM cannot be located, and we want it as a submodule")
  endif()
  if(CMAKE_CROSSCOMPILING)
    set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
  endif()
  find_package(IFEMDeps) # To get FOUND variables in current context
  add_subdirectory(${IFEM_PATH} IFEM)
  include(${CMAKE_BINARY_DIR}/IFEM/IFEMFlags.cmake)
  set(IFEM_LIBRARIES IFEM)
else()
  # Build is in-tree
  message(STATUS "Using in-tree IFEM")
  set(IFEM_INTREE_BUILD ON)
  set(IFEM_INCLUDES ${IFEM_H_PATH}
                    ${IFEM_PATH}/src/ASM
                    ${IFEM_PATH}/src/Eig
                    ${IFEM_PATH}/src/LinAlg
                    ${IFEM_PATH}/src/SIM
                    ${IFEM_PATH}/src/Utility
                    ${IFEM_PATH}/3rdparty/gtest/include)

  find_library(IFEM_LIBRARIES
    NAMES IFEM
    PATHS ${IFEM_PATH}/${IFEM_BUILD_TYPE}/lib
          ${IFEM_PATH}/lib
    NO_DEFAULT_PATH)
  if(NOT IFEM_LIBRARIES)
    message(WARNING "Could not find the in-tree IFEM library, "
                    "we assume it will be built into a build-type dir")
  endif()
  if(NOT BUILD_SHARED_LIBS)
    find_package(IFEMDeps REQUIRED)
    set(IFEM_LIBRARIES ${IFEM_LIBRARIES} ${IFEM_DEPLIBS})
  endif()
  add_library(IFEMAppCommon STATIC IMPORTED)
  set_target_properties(IFEMAppCommon PROPERTIES IMPORTED_LOCATION
                        ${IFEM_PATH}/${IFEM_BUILD_TYPE}/lib/libIFEMAppCommon${CMAKE_STATIC_LIBRARY_SUFFIX})
  if(LRSpline_FOUND)
    list(APPEND IFEM_INCLUDES ${IFEM_PATH}/src/ASM/LR)
  endif()
  if(NOT IFEM_USE_SYSTEM_TINYXML)
    set(IFEM_INCLUDES ${IFEM_INCLUDES}
                      ${IFEM_PATH}/3rdparty/tinyxml)
  endif()
  set(IFEM_LIBRARIES ${IFEM_LIBRARIES} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif()

list(APPEND CMAKE_MODULE_PATH ${IFEM_PATH}/cmake/Scripts)

if(NOT IFEM_AS_SUBMODULE)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(IFEM DEFAULT_MSG
                                    IFEM_LIBRARIES)

  # Export version information
  if(IFEM_INCLUDES)
    execute_process(COMMAND cat "${IFEM_H_PATH}/IFEM.h" OUTPUT_VARIABLE IFEM_HEADER)
    string(REGEX MATCH "IFEM_VERSION_MAJOR ([0-9]+)" IFEM_VERSION_MAJOR ${IFEM_HEADER})
    string(REGEX REPLACE "IFEM_VERSION_MAJOR ([0-9]+)" "\\1" IFEM_VERSION_MAJOR "${IFEM_VERSION_MAJOR}")
    string(REGEX MATCH "IFEM_VERSION_MINOR ([0-9]+)" IFEM_VERSION_MINOR ${IFEM_HEADER})
    string(REGEX REPLACE "IFEM_VERSION_MINOR ([0-9]+)" "\\1" IFEM_VERSION_MINOR "${IFEM_VERSION_MINOR}")
    string(REGEX MATCH "IFEM_VERSION_PATCH ([0-9]+)" IFEM_VERSION_PATCH ${IFEM_HEADER})
    string(REGEX REPLACE "IFEM_VERSION_PATCH ([0-9]+)" "\\1" IFEM_VERSION_PATCH "${IFEM_VERSION_PATCH}")
    set(IFEM_VERSION "${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR}.${IFEM_VERSION_PATCH}")
    set(IFEM_ABI_VERSION ${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR})
  endif()

  set(IFEM_LIBRARIES ${IFEM_LIBRARIES}
                     ${IFEM_DEPLIBS})
  set(IFEM_INCLUDES ${IFEM_INCLUDES} ${IFEM_DEPINCLUDES})
endif()

if(IFEM_TEST_MEMCHECK)
  include(CTest)
endif()

# Needed as we have templates using these flags
enable_language(CXX)
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DUSE_MKL -mkl=sequential")
  set(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DUSE_MKL -mkl=sequential")
  list(APPEND IFEM_DEFINITIONS -DUSE_MKL=1)
else()
  set(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DUSE_CBLAS")
  set(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DUSE_CBLAS")
  list(APPEND IFEM_DEFINITIONS -DUSE_CBLAS=1)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DINDEX_CHECK=2")
if(VERBOSE_DEBUG GREATER 0)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DINT_DEBUG=${VERBOSE_DEBUG}")
endif()
set(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DReal=double")
list(APPEND IFEM_DEFINITIONS -DReal=double)

set(IFEM_CONFIGURED 1)
if(NOT IFEM_TESTING_INCLUDED)
  include(../Scripts/IFEMTesting)
endif()

include(../Scripts/IFEMDoxy)
