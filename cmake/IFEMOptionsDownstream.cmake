# Common options for IFEM applications

option(IFEM_AS_SUBMODULE       "Compile IFEM as a submodule of apps?"       OFF)
option(IFEM_TEST_MEMCHECK      "Run tests through valgrind?"                OFF)
option(IFEM_BUILD_TESTING      "Build testing by default?"                  OFF)
option(IFEM_SERIAL_TESTS_IN_PARALLEL "Run serial tests in parallel builds?" ON)
option(IFEM_INSTALL_DOXY       "Install documentation?"                     ON)
option(IFEM_USE_NATIVE         "Enable native tuning in release builds?"    ON)
