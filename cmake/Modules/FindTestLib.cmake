find_package(GTest QUIET)
if(NOT TestLib_FOUND)
  if(GTest_FOUND)
    set(TESTLIB_LIBRARY GTest::GTest)
  else()
    find_package(Threads REQUIRED)
    if(IFEM_INTREE_BUILD OR IFEM_LIBRARY_BUILD OR IFEM_AS_SUBMODULE)
       add_subdirectory(${IFEM_PATH}/3rdparty/gtest gtest EXCLUDE_FROM_ALL)
       add_library(GTest::GTest UNKNOWN IMPORTED)
       set_target_properties(GTest::GTest PROPERTIES
                             IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/gtest/libgtest.a
                             INTERFACE_LINK_LIBRARIES Threads::Threads
                             INTERFACE_INCLUDE_DIRECTORIES ${IFEM_PATH}/3rdparty/gtest/include)
       add_dependencies(GTest::GTest gtest)
     else()
      include(DownloadGTest)
    endif()
    get_target_property(TESTLIB_LIBRARY GTest::GTest IMPORTED_LOCATION)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TestLib DEFAULT_MSG
                                  TESTLIB_LIBRARY)
if(TestLib_FOUND)
  set(GTest_FOUND CACHE BOOL "" FORCE)
endif()
