project(IFEMAppCommon)

cmake_minimum_required(VERSION 2.6)

include_directories(${PROJECT_SOURCE_DIR})

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

if(IFEM_AS_SUBMODULE)
  set(IFEM_LIBRARIES IFEM ${IFEM_DEPLIBS})
endif()

# Common application sources
add_library(IFEMAppCommon STATIC MeshUtils.C
                                 Spalding.C
                                 SAWallLaw.C
                                 SIMNodalConstraint.C
                                 SIMSemi3D.C
                                 StabilizationUtils.C
                                 TimeIntUtils.C
                                 WeakOperators.C)
install(TARGETS IFEMAppCommon DESTINATION ${CMAKE_INSTALL_LIBDIR})
file(GLOB AppCommon_HEADERS *.h)
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
set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
