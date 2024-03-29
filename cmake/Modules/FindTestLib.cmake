# Shipped cmake modules for gtest relies on policy to avoid warnings
cmake_policy(SET CMP0012 NEW)

find_package(GTest QUIET)
if(NOT TestLib_FOUND)
  if(GTest_FOUND)
    set(TESTLIB_LIBRARY GTest::GTest)
  else()
    find_package(Threads REQUIRED)
    include(DownloadGTest)
    get_target_property(TESTLIB_LIBRARY GTest::GTest IMPORTED_LOCATION)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TestLib DEFAULT_MSG
                                  TESTLIB_LIBRARY)
if(TestLib_FOUND)
  set(GTest_FOUND CACHE BOOL "" FORCE)
endif()
