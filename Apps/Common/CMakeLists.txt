cmake_minimum_required(VERSION 3.5)

project(IFEMAppCommon)

include_directories(${PROJECT_SOURCE_DIR})

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(IFEM_AS_SUBMODULE)
  set(IFEM_LIBRARIES IFEM ${IFEM_DEPLIBS})
endif()

# Common application sources
file(GLOB AppCommon_SOURCES *.C)
file(GLOB AppCommon_HEADERS *.h)
add_library(IFEMAppCommon STATIC ${AppCommon_SOURCES})
install(TARGETS IFEMAppCommon DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${AppCommon_HEADERS} DESTINATION include/IFEM)
install(FILES scripts/regtest.sh.in
              scripts/clang-check-test.sh.in
              scripts/cppcheck-test.sh.in
              scripts/iotest.sh.in
        DESTINATION lib/IFEM/scripts
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE)

# Unit tests
file(GLOB APPCOMMON_TEST_SRCS ${PROJECT_SOURCE_DIR}/Test/*.C)
if(LRSpline_FOUND)
  list(APPEND APPCOMMON_TEST_SRCS ${PROJECT_SOURCE_DIR}/Test/LR/TestMultiPatchLRRefine.C)
endif()

IFEM_add_test_app("${APPCOMMON_TEST_SRCS}"
                  ${PROJECT_SOURCE_DIR}/Test
                  AppCommon 0
                  IFEMAppCommon ${IFEM_LIBRARIES})

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
list(APPEND CHECK_SOURCES ${AppCommon_SOURCES})
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
