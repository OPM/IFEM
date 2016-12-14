project(IFEMAppCommon)

cmake_minimum_required(VERSION 2.6)

include_directories(${PROJECT_SOURCE_DIR})

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

if(IFEM_AS_SUBMODULE)
  set(IFEM_LIBRARIES IFEM ${IFEM_DEPLIBS})
endif()

# Common application sources
file(GLOB AppCommon_SOURCES *.C)
file(GLOB AppCommon_HEADERS *.h)
add_library(IFEMAppCommon STATIC ${AppCommon_SOURCES})
install(TARGETS IFEMAppCommon DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${AppCommon_HEADERS} DESTINATION include/IFEM)
install(FILES scripts/regtest.sh.in DESTINATION lib/IFEM/scripts
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  AppCommon
                  IFEMAppCommon ${IFEM_LIBRARIES})
if(MPI_FOUND)
  set(TEST_SRCS_MPI ${PROJECT_SOURCE_DIR}/Test/MPI/TestMultiPatchModelGenerator.C)
  add_executable(AppCommon-MPI-test EXCLUDE_FROM_ALL
                 ${IFEM_PATH}/src/IFEM-test.C ${TEST_SRCS_MPI})
  target_link_libraries(AppCommon-MPI-test IFEMAppCommon ${IFEM_LIBRARIES} ${IFEM_DEPLIBS} gtest)
  gtest_add_tests($<TARGET_FILE:AppCommon-MPI-test> ${IFEM_PATH}/Apps/Common TEST_SRCS_MPI
                  ${MPIEXEC} -np 4)
  list(APPEND TEST_APPS AppCommon-MPI-test)
endif()

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
