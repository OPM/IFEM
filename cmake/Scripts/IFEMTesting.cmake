macro(IFEM_add_test_app path workdir name parallel)
  if("${path}" MATCHES "\\*")
    file(GLOB TEST_SRCS ${path})
  else()
    set(TEST_SRCS ${path})
  endif()
  if(IFEM_BUILD_TESTING)
    set(EXCL_ALL)
  else()
    set(EXCL_ALL EXCLUDE_FROM_ALL)
  endif()
  add_executable(${name}-test ${EXCL_ALL} ${IFEM_PATH}/src/IFEM-test.C ${TEST_SRCS})
  if(${parallel} GREATER 0)
    set_property(TARGET ${name}-test PROPERTY
                 CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${parallel}')
  endif()
  include(GoogleTest)
  gtest_discover_tests(${name}-test
                       WORKING_DIRECTORY ${workdir}
                       NO_PRETTY_VALUES)
  list(APPEND TEST_APPS ${name}-test)
  target_link_libraries(${name}-test GTest::GTest ${ARGN})
endmacro()

macro(IFEM_add_unittests IFEM_PATH)
  file(GLOB TEST_SOURCES ${IFEM_PATH}/src/Utility/Test/*.C;${IFEM_PATH}/src/ASM/Test/*.C;${IFEM_PATH}/src/LinAlg/Test/*.C;${IFEM_PATH}/src/SIM/Test/*.C;${IFEM_PATH}/src/Eig/Test/*.C)
  if(NOT PETSC_FOUND)
    list(REMOVE_ITEM TEST_SOURCES ${IFEM_PATH}/src/LinAlg/Test/TestPETScMatrix.C)
  endif()

  if(NOT ISTL_FOUND)
    list(REMOVE_ITEM TEST_SOURCES ${IFEM_PATH}/src/LinAlg/Test/TestISTLMatrix.C)
  endif()

  if(NOT HDF5_FOUND)
    list(REMOVE_ITEM TEST_SOURCES ${IFEM_PATH}/src/Utility/Test/TestFieldFunctions.C)
    list(REMOVE_ITEM TEST_SOURCES ${IFEM_PATH}/src/Utility/Test/TestFieldFunctionsLR.C)
  endif()

  if(LRSPLINE_FOUND OR LRSpline_FOUND)
    file(GLOB LR_TEST_SRCS ${IFEM_PATH}/src/ASM/LR/Test/*.C)
    list(APPEND TEST_SOURCES ${LR_TEST_SRCS})
  else()
    list(REMOVE_ITEM TEST_SOURCES ${IFEM_PATH}/src/Utility/Test/TestFieldFunctionsLR.C)
  endif()

  IFEM_add_test_app("${TEST_SOURCES}"
                    ${IFEM_PATH}
                    IFEM 0
                    ${IFEM_LIBRARIES} ${IFEM_DEPLIBS})

  IFEM_add_test_app("${IFEM_PATH}/src/LinAlg/Test/NoBlas/TestMatrixFallback.C"
                    ${IFEM_PATH}
                    IFEM_noblas 0 IFEM)

  # Parallel unit tests. These all run with 4 processes.
  if(MPI_FOUND)
    set(TEST_SRCS_MPI ${IFEM_PATH}/src/ASM/Test/MPI/TestDomainDecomposition.C)
    if(PETSC_FOUND)
      list(APPEND TEST_SRCS_MPI ${IFEM_PATH}/src/LinAlg/Test/MPI/TestPETScMatrix.C)
    endif()
    if(ISTL_FOUND)
      list(APPEND TEST_SRCS_MPI ${IFEM_PATH}/src/LinAlg/Test/MPI/TestISTLMatrix.C)
    endif()
    IFEM_add_test_app("${TEST_SRCS_MPI}" ${IFEM_PATH} IFEM-MPI 4 ${IFEM_LIBRARIES} ${IFEM_DEPLIBS})
    list(APPEND TEST_APPS IFEM-MPI-test)
  endif()
endmacro()

function(IFEM_add_test name binary)
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_EXTRA)
    set(test-name "${binary}+${IFEM_TEST_EXTRA}+${name}")
  else()
    set(test-name "${binary}+${name}")
  endif()
  if(ARGN)
    set(ARGN "-m ${ARGN}")
  endif()
  if(IFEM_TEST_MEMCHECK)
    set(VALGRIND "-v ${MEMCHECK_COMMAND}")
  endif()
  add_test("${test-name}" regtest.sh ${VALGRIND} -a ${EXECUTABLE_OUTPUT_PATH}/${binary} -r ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${name} ${ARGN})
  set_tests_properties(${test-name} PROPERTIES TIMEOUT 5000)
endfunction()

function(IFEM_add_restart_test name binary rlevel)
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_EXTRA)
    set(test-name "restart+${binary}+${IFEM_TEST_EXTRA}+${name}")
  else()
    set(test-name "restart+${binary}+${name}")
  endif()
  if(ARGN)
    set(ARGN "-m ${ARGN}")
  endif()
  if(IFEM_TEST_MEMCHECK)
    set(VALGRIND "-v ${MEMCHECK_COMMAND}")
  endif()
  add_test("${test-name}" regtest.sh ${VALGRIND} -a ${EXECUTABLE_OUTPUT_PATH}/${binary} -r ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${name} -s ${rlevel} ${ARGN})
endfunction()

function(IFEM_add_vtf_test name binary)
  if(NOT VTFWRITER_FOUND OR NOT VTFLS_COMMAND)
    return()
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_EXTRA)
    set(test-name "vtf+${binary}+${IFEM_TEST_EXTRA}+${name}")
  else()
    set(test-name "vtf+${binary}+${name}")
  endif()
  add_test("${test-name}" iotest.sh ${EXECUTABLE_OUTPUT_PATH}/${binary} ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${name} vtf)
endfunction()

function(IFEM_add_hdf5_test name binary)
  if(NOT HDF5_FOUND OR NOT H5LS_COMMAND)
    return()
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_EXTRA)
    set(test-name "hdf5+${binary}+${IFEM_TEST_EXTRA}+${name}")
  else()
    set(test-name "hdf5+${binary}+${name}")
  endif()
  add_test("${test-name}" iotest.sh ${EXECUTABLE_OUTPUT_PATH}/${binary} ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${name} hdf5 ${ARGN})
endfunction()

macro(add_check_target)
  add_custom_target(check ${CMAKE_CTEST_COMMAND} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  add_custom_command(TARGET check PRE_BUILD COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/failed.log)
  if(IFEM_AS_SUBMODULE OR IFEM_LIBRARY_BUILD)
    ifem_add_unittests(${IFEM_PATH})
  endif()
  add_dependencies(check ${TEST_APPS})
  add_custom_target(testapps DEPENDS ${TEST_APPS})

  if(NOT TARGET check-commits)
    add_custom_target(check-commits
                      COMMAND ${CMAKE_COMMAND}
                              -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                              -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
                              -P ${IFEM_CHECKCOMMITS_SCRIPT})
  endif()

  foreach(dep ${IFEM_INCLUDE_DIRS} ${IFEM_INCLUDES})
    list(APPEND IPATHS -I ${dep})
  endforeach()

  if(NOT CMAKE_DISABLE_FIND_PACKAGE_ClangCheck)
    find_program(CLANGCHECK_COMMAND clang-check)
  endif()
  if(NOT CMAKE_DISABLE_FIND_PACKAGE_CppCheck)
    find_program(CPPCHECK_COMMAND cppcheck)
  endif()
  foreach(src ${CHECK_SOURCES})
    get_filename_component(name ${src} NAME)
    get_filename_component(EXT ${src} EXT)
    if(EXT STREQUAL .C)
      if(NOT IS_ABSOLUTE ${src})
        set(src ${PROJECT_SOURCE_DIR}/${src})
      endif()
      if(CLANGCHECK_COMMAND AND CMAKE_EXPORT_COMPILE_COMMANDS)
        add_test(NAME clang-check+${name}
                 COMMAND clang-check-test.sh ${CLANGCHECK_COMMAND} ${src}
                 CONFIGURATIONS analyze clang-check)
      endif()
      if(CPPCHECK_COMMAND)
        add_test(NAME cppcheck+${name}
                 COMMAND cppcheck-test.sh ${CPPCHECK_COMMAND} ${src} ${IPATHS}
                 CONFIGURATIONS analyze cppcheck)
      endif()
    endif()
  endforeach()
endmacro()

# For regression tests
if(IFEM_TEST_MEMCHECK)
  find_program(MEMCHECK_COMMAND valgrind)
  if(NOT MEMCHECK_COMMAND)
    message(FATAL_ERROR "Could not locate valgrind and IFEM_TEST_MEMCHECK is requested.")
  endif()
endif()

# Used for unit tests
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=yes")
include(CTest)

find_package(TestLib REQUIRED)

find_program(VTFLS_COMMAND vtfls)
find_program(H5LS_COMMAND h5ls)

# Generate regtest script with correct paths
configure_file(${IFEM_REGTEST_SCRIPT} regtest.sh @ONLY)
configure_file(${IFEM_CLANG_CHECK_TEST_SCRIPT} clang-check-test.sh @ONLY)
configure_file(${IFEM_CPPCHECK_TEST_SCRIPT} cppcheck-test.sh @ONLY)
configure_file(${IFEM_IOTEST_SCRIPT} iotest.sh @ONLY)

set(IFEM_TESTING_INCLUDED 1 CACHE BOOL "" FORCE)
