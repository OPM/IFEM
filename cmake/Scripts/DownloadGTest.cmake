########################### GTEST
# Enable ExternalProject CMake module
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/third_party)

# Add gtest
# http://stackoverflow.com/questions/9689183/cmake-googletest
externalproject_add(
    googletest
    URL http://github.com/google/googletest/archive/release-1.8.0.tar.gz
    BUILD_BYPRODUCTS third_party/src/googletest-build/googlemock/gtest/libgtest.a
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

# Specify include dir
externalproject_get_property(googletest source_dir)
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)

# Library
externalproject_get_property(googletest binary_dir)
set(GTEST_LIBRARIES ${binary_dir}/googlemock/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
add_library(gtest UNKNOWN IMPORTED)
set_property(TARGET gtest PROPERTY IMPORTED_LOCATION
             ${GTEST_LIBRARIES})
set_property(TARGET googletest PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(gtest googletest)
