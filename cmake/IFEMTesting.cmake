# Add a catch2 unit test application
#   Registers tests as ctest entries using test discovery mechanism.
# Single-valued parameters:
#   NAME      - Name of test application
#     The target added will be ${NAME]-test.
#   PARALLEL  - Number of MPI procs to use.
#     If given, tests are only added if MPI is enabled.
#     If not given, tests are not added if MPI is enabled
#     and IFEM_SERIAL_TESTS_IN_PARALLEL is OFF.
#   WORKDIR   - Directory tests will be executed in
# Multi-valued parameters:
#   LIBRARIES - Libraries test application will be linked to
#   SOURCES   - Source files for tests
function(ifem_add_test_app)
  set(oneValueArgs NAME PARALLEL WORKDIR)
  set(multiValueArgs LIBRARIES SOURCES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(PARAM_PARALLEL)
    if(NOT MPI_FOUND)
      return()
    endif()
  else()
    if(MPI_FOUND AND NOT IFEM_SERIAL_TESTS_IN_PARALLEL)
      return()
    endif()
  endif()

  foreach(source ${PARAM_SOURCES})
    if(source MATCHES "\\*")
      file(GLOB TEST_SRCS ${source})
      list(APPEND test_sources ${TEST_SRCS})
    else()
      list(APPEND test_sources ${source})
    endif()
  endforeach()
  if(IFEM_BUILD_TESTING)
    set(EXCL_ALL)
  else()
    set(EXCL_ALL EXCLUDE_FROM_ALL)
  endif()
  add_executable(${PARAM_NAME}-test ${EXCL_ALL} ${IFEM_PATH}/src/IFEM-test.C ${test_sources})
  set_target_properties(
    ${PARAM_NAME}-test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/bin
  )
  if(PARAM_PARALLEL GREATER 0)
    set_property(
      TARGET ${PARAM_NAME}-test
      PROPERTY
        CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PARAM_PARALLEL}'
    )
    catch_discover_tests(
      ${PARAM_NAME}-test
      WORKING_DIRECTORY
        ${PARAM_WORKDIR}
      PROPERTIES
        PROCESSORS ${PARAM_PARALLEL}
        LABELS unit
      NO_PRETTY_VALUES
    )
  else()
    catch_discover_tests(
      ${PARAM_NAME}-test
      WORKING_DIRECTORY
        ${PARAM_WORKDIR}
      PROPERTIES
        LABELS unit
      NO_PRETTY_VALUES
    )
  endif()
  list(APPEND TEST_APPS ${PARAM_NAME}-test)
  find_package(Catch2 REQUIRED)
  target_link_libraries(${PARAM_NAME}-test PRIVATE Catch2::Catch2 ${PARAM_LIBRARIES})
  target_compile_definitions(${PARAM_NAME}-test PRIVATE CATCH2_VERSION_MAJOR=${Catch2_VERSION_MAJOR})
  ifem_add_sca_tests(TARGET ${PARAM_NAME}-test)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
endfunction()

# Add regression tests
#   Handles enabling/disabling tests according to MPI configuration.
# Single-valued parameters:
#   PARALLEL   - Number of MPI procs to use.
#     If given, tests are only added if MPI is enabled.
#     If not given, tests are not added if MPI is enabled
#     and IFEM_SERIAL_TESTS_IN_PARALLEL is OFF.
#   TARGET     - Target for application to add regression test for
# Multi-valued parameters:
#   CONDITIONS - Conditions for enabling the tests
#   TEST_FILES - Name of .reg files for tests
function(ifem_add_regression_test)
  set(oneValueArgs PARALLEL TARGET)
  set(multiValueArgs CONDITIONS TEST_FILES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(PARAM_CONDITIONS)
    foreach(condition ${PARAM_CONDITIONS})
      if(NOT (${condition}))
        return()
      endif()
    endforeach()
  endif()
  if(PARAM_PARALLEL)
    if(NOT MPI_FOUND)
      return()
    endif()
    set(PARALLEL "-m ${PARAM_PARALLEL}")
  else()
    if(MPI_FOUND AND NOT IFEM_SERIAL_TESTS_IN_PARALLEL)
      return()
    endif()
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_MEMCHECK)
    set(VALGRIND -v ${MEMCHECK_COMMAND})
  endif()
  foreach(regfile ${PARAM_TEST_FILES})
    set(test_name "${PARAM_TARGET}+${regfile}")
    add_test(NAME ${test_name}
            COMMAND
            ifem-reg-test
            -a $<TARGET_FILE:${PARAM_TARGET}>
            -r ${PROJECT_SOURCE_DIR}/Test/${regfile}
            -d ${CMAKE_BINARY_DIR}
            ${PARALLEL}
            ${VALGRIND}
    )
    set_tests_properties(
      ${test_name}
      PROPERTIES
        TIMEOUT
          10000
        LABELS
          regression
    )
    if(PARAM_PARALLEL)
      set_tests_properties(${test_name} PROPERTIES PROCESSORS ${PARAM_PARALLEL})
    endif()
  endforeach()
endfunction()

# Add restart tests
#   Handles enabling/disabling tests according to MPI, HDF5 and cereal configuration.
# Single-valued parameters:
#   PARALLEL      - Number of MPI procs to use.
#     If given, tests are only added if MPI is enabled.
#     If not given, tests are not added if MPI is enabled
#     and IFEM_SERIAL_TESTS_IN_PARALLEL is OFF.
#   RESTART_LEVEL - Time level to restart at
#   TARGET        - Target for application to add regression test for
# Multi-valued parameters:
#   CONDITIONS    - Conditions for enabling the tests
#   TEST_FILES    - Name of .reg files for tests
function(ifem_add_restart_test)
  if(NOT cereal_FOUND OR NOT HDF5_FOUND)
    return()
  endif()
  set(oneValueArgs TARGET PARALLEL RESTART_LEVEL)
  set(multiValueArgs CONDITIONS TEST_FILES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(PARAM_CONDITIONS)
    foreach(condition ${PARAM_CONDITIONS})
      if(NOT (${condition}))
        return()
      endif()
    endforeach()
  endif()
  if(PARAM_PARALLEL)
    if(NOT MPI_FOUND)
      return()
    endif()
    set(PARALLEL "-m ${PARAM_PARALLEL}")
  else()
    if(MPI_FOUND AND NOT IFEM_SERIAL_TESTS_IN_PARALLEL)
      return()
    endif()
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_MEMCHECK)
    set(VALGRIND -v ${MEMCHECK_COMMAND})
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  foreach(regfile ${PARAM_TEST_FILES})
    set(test_name "restart+${PARAM_TARGET}+${regfile}")
    add_test(NAME ${test_name}
            COMMAND
            ifem-reg-test
            -a $<TARGET_FILE:${PARAM_TARGET}>
            -r ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${regfile}
            -s ${PARAM_RESTART_LEVEL}
            -d ${CMAKE_BINARY_DIR}
            ${VALGRIND}
            ${PARALLEL}
    )
    set_tests_properties(
      ${test_name}
      PROPERTIES
        TIMEOUT
          10000
        LABELS
          restart
    )
    if(PARAM_PARALLEL)
      set_tests_properties(${test_name} PROPERTIES PROCESSORS ${PARAM_PARALLEL})
    endif()
  endforeach()
endfunction()

# Add VTF tests
#   Handles enabling/disabling tests according to VTF configuration.
# Single-valued parameters:
#   TARGET        - Target for application to add regression test for
# Multi-valued parameters:
#   CONDITIONS    - Conditions for enabling the tests
#   TEST_FILES    - Name of .vreg files for tests
function(ifem_add_vtf_test)
  if(NOT TARGET vtfls OR NOT VTFWriter_FOUND)
    return()
  endif()
  if(MPI_FOUND AND NOT IFEM_SERIAL_TESTS_IN_PARALLEL)
    return()
  endif()
  set(oneValueArgs TARGET)
  set(multiValueArgs CONDITIONS TEST_FILES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(PARAM_CONDITIONS)
    foreach(condition ${PARAM_CONDITIONS})
      if(NOT (${condition}))
        return()
      endif()
    endforeach()
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_MEMCHECK)
    set(VALGRIND -v ${MEMCHECK_COMMAND})
  endif()
  foreach(regfile ${PARAM_TEST_FILES})
    set(test_name "vtf+${PARAM_TARGET}+${regfile}")
    add_test(NAME ${test_name}
             COMMAND
             ifem-io-test
             -a $<TARGET_FILE:${PARAM_TARGET}>
             -r ${PROJECT_SOURCE_DIR}/Test/${regfile}
             -t vtf
             -p $<TARGET_FILE:vtfls>
             -d ${CMAKE_BINARY_DIR}
             ${VALGRIND}
    )
    set_tests_properties(
      ${test_name}
      PROPERTIES
        TIMEOUT
          10000
        LABELS
          "io;vtf"
    )
  endforeach()
endfunction()

# Add HDF5 tests
#   Handles enabling/disabling tests according to HDF5 and MPI configuration.
# Single-valued parameters:
#   PARALLEL      - Number of MPI procs to use.
#     If given, tests are only added if MPI is enabled.
#     If not given, tests are not added if MPI is enabled
#     and IFEM_SERIAL_TESTS_IN_PARALLEL is OFF.
#   TARGET        - Target for application to add regression test for
# Multi-valued parameters:
#   CONDITIONS    - Conditions for enabling the tests
#   TEST_FILES    - Name of .hreg files for tests
function(ifem_add_hdf5_test)
  if(NOT TARGET h5ls OR NOT HDF5_FOUND)
    return()
  endif()
  set(oneValueArgs TARGET PARALLEL)
  set(multiValueArgs TEST_FILES CONDITIONS)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(PARAM_PARALLEL)
    if(NOT MPI_FOUND)
      return()
    endif()
    set(PARAM_PARALLEL "-m ${PARAM_PARALLEL}")
  else()
    if(MPI_FOUND AND NOT IFEM_SERIAL_TESTS_IN_PARALLEL)
      return()
    endif()
  endif()
  if(PARAM_CONDITIONS)
    foreach(condition ${PARAM_CONDITIONS})
      if(NOT (${condition}))
        return()
      endif()
    endforeach()
  endif()
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_MEMCHECK)
    set(VALGRIND -v ${MEMCHECK_COMMAND})
  endif()
  foreach(regfile ${PARAM_TEST_FILES})
    set(test_name "hdf5+${PARAM_TARGET}+${regfile}")
    add_test(NAME ${test_name}
             COMMAND
             ifem-io-test
             -a $<TARGET_FILE:${PARAM_TARGET}>
             -r ${PROJECT_SOURCE_DIR}/Test/${regfile}
             -t hdf5
             -p $<TARGET_FILE:h5ls>
             -d ${CMAKE_BINARY_DIR}
             ${PARAM_PARALLEL}
             ${VALGRIND}
    )
    set_tests_properties(
      ${test_name}
      PROPERTIES
        TIMEOUT
          10000
        LABELS
          "io;hdf5"
    )
    if(PARAM_PARALLEL)
      set_tests_properties(${test_name} PROPERTIES PROCESSORS ${PARAM_PARALLEL})
    endif()
  endforeach()
endfunction()

# Add check and testapps convenience targets
macro(ifem_add_check_target)
  add_custom_target(check ${CMAKE_CTEST_COMMAND} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  add_custom_command(TARGET check PRE_BUILD COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/failed.log)
  add_dependencies(check ${TEST_APPS})
  add_custom_target(testapps DEPENDS ${TEST_APPS})
endmacro()

# Add static analysis tests
#  Adds clang-check tests to ctest according to clang-check availability.
# Multi-valued parameters:
#   TARGETS - Targets for source files to add tests for
function(ifem_add_sca_tests)
  if(NOT CMAKE_DISABLE_FIND_PACKAGE_ClangCheck)
    find_program(CLANGCHECK_COMMAND clang-check)
  endif()
  set(multiValueArgs TARGET)
  cmake_parse_arguments(PARAM "" "" "${multiValueArgs}" ${ARGN})

  foreach(target ${PARAM_TARGET})
    get_target_property(SOURCES ${target} SOURCES)

    foreach(src ${SOURCES})
      get_filename_component(name ${src} NAME)
      get_filename_component(EXT ${src} EXT)
      if(EXT STREQUAL .C)
        if(NOT IS_ABSOLUTE ${src})
          set(src ${PROJECT_SOURCE_DIR}/${src})
        endif()
        if(CLANGCHECK_COMMAND AND CMAKE_EXPORT_COMPILE_COMMANDS)
          if(NOT TEST clang-check+${name})
            add_test(
              NAME
                clang-check+${name}
              COMMAND
                ifem-clang-check-test
                ${CLANGCHECK_COMMAND}
                ${src}
                ${PROJECT_BINARY_DIR}
              CONFIGURATIONS
                analyze
                clang-check
            )
          endif()
        endif()
      endif()
    endforeach()
  endforeach()
endfunction()

if(IFEM_TEST_MEMCHECK)
  find_program(MEMCHECK_COMMAND valgrind)
  if(NOT MEMCHECK_COMMAND)
    message(FATAL_ERROR
      "Could not locate valgrind and IFEM_TEST_MEMCHECK is requested."
    )
  endif()
endif()

find_package(Catch2 REQUIRED)
include(Catch)

include(CTest)

find_program(VTFLS_COMMAND vtfls)
if(VTFLS_COMMAND)
  add_executable(vtfls IMPORTED GLOBAL)
  set_target_properties(vtfls PROPERTIES IMPORTED_LOCATION ${VTFLS_COMMAND})
endif()

find_program(H5LS_COMMAND h5ls)
if(H5LS_COMMAND)
  add_executable(h5ls IMPORTED GLOBAL)
  set_target_properties(h5ls PROPERTIES IMPORTED_LOCATION ${H5LS_COMMAND})
endif()

add_executable(ifem-reg-test IMPORTED GLOBAL)
set_target_properties(ifem-reg-test PROPERTIES
                       IMPORTED_LOCATION
                         ${IFEM_REGTEST_PATH}
)
add_executable(ifem-clang-check-test IMPORTED GLOBAL)
set_target_properties(ifem-clang-check-test PROPERTIES
                      IMPORTED_LOCATION
                        ${IFEM_CLANG_CHECK_PATH}
)
add_executable(ifem-io-test IMPORTED GLOBAL)
set_target_properties(ifem-io-test PROPERTIES
                      IMPORTED_LOCATION
                        ${IFEM_IOTEST_PATH}
)

set(IFEM_TESTING_INCLUDED 1 CACHE BOOL "" FORCE)
