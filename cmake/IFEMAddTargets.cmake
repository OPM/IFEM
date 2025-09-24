# Function for setting up a static library
#   Adds a static library, setting up sources, filesets and associated tests
# Single-valued parameters:
#   NAME      - Name of library
#     The target added will be ${NAME}
# Multi-valued parameters:
#   BASE_DIRS - Base directories for the headers file set
#   HEADERS   - Header files associated with library
#   LIBRARIES - Targets to link the library to
#   SOURCES   - Source files for the library
function(ifem_add_library)
  set(oneValueArgs NAME)
  set(multiValueArgs BASE_DIRS HEADERS LIBRARIES SOURCES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  add_library(${PARAM_NAME} STATIC)
  target_sources(${PARAM_NAME} PRIVATE ${PARAM_SOURCES})
  if(PARAM_HEADERS)
    if(CMAKE_VERSION GREATER_EQUAL 3.23)
      target_sources(
        ${PARAM_NAME} PUBLIC
        FILE_SET
          HEADERS
        BASE_DIRS
          ${PROJECT_SOURCE_DIR}
          ${PARAM_BASE_DIRS}
        FILES
          ${PARAM_HEADERS})
    else()
      foreach(dir ${PROJECT_SOURCE_DIR} ${PARAM_BASE_DIRS})
        target_include_directories(${PARAM_NAME} PUBLIC $<BUILD_INTERFACE:${dir}>)
      endforeach()
    endif()
  endif()
  target_link_libraries(${PARAM_NAME} PUBLIC ${PARAM_LIBRARIES})
  set_target_properties(
    ${PARAM_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/lib
  )
  ifem_add_sca_tests(TARGET ${PARAM_NAME})
endfunction()

# Function for setting up an application
#   Adds an executable, setting up sources, filesets and associated tests
# Single-valued parameters:
#   NAME      - Name of application
#     The target added will be ${NAME}
# Multi-valued parameters:
#   HEADERS   - Header files associated with library
#   LIBRARIES - Targets to link the library to
#   SOURCES   - Source files for the library
function(ifem_add_application)
  set(oneValueArgs NAME)
  set(multiValueArgs HEADERS LIBRARIES SOURCES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  add_executable(${PARAM_NAME})
  target_sources(${PARAM_NAME} PRIVATE ${PARAM_SOURCES})
  if(PARAM_HEADERS)
    if(CMAKE_VERSION GREATER_EQUAL 3.23)
      target_sources(${PARAM_NAME} PUBLIC FILE_SET HEADERS FILES ${PARAM_HEADERS})
    endif()
  endif()
  target_link_libraries(${PARAM_NAME} PUBLIC ${PARAM_LIBRARIES})
  set_target_properties(
    ${PARAM_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/bin
  )
  ifem_add_sca_tests(TARGET ${PARAM_NAME})
  list(APPEND TEST_APPS ${PARAM_NAME})
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
endfunction()

# Function for adding a dependency dir
#   Add subdirectory from application, handling the IFEM- prefix ambiguity.
# No-valued parameters:
#   OPTIONAL - True if dependency is optional
# Single-valued parameters:
#   APP      - Name of application to add directory from
#   DEPTH    - How many directory entries to remove for locating application
#   TARGET   - Target we expect to import from subdirectory
# Multi-valued parameters:
#   DIR      - Possible names for subdirectory to add
function(ifem_add_dependency_dir)
  set(options OPTIONAL)
  set(oneValueArgs APP DEPTH TARGET)
  set(multiValueArgs DIR)
  cmake_parse_arguments(PARAM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT PARAM_DEPTH)
    set(PARAM_DEPTH 1)
  endif()

  string(REPEAT "/.." ${PARAM_DEPTH} DOTS)
  set(TOP_DIR ${PROJECT_SOURCE_DIR}${DOTS})

  if(NOT TARGET ${PARAM_TARGET})
    if (EXISTS ${TOP_DIR}/IFEM-${PARAM_APP})
      set(APP_DIR ${TOP_DIR}/IFEM-${PARAM_APP})
    elseif(EXISTS ${TOP_DIR}/${PARAM_APP})
      set(APP_DIR ${TOP_DIR}/${PARAM_APP})
    else()
      if(NOT PARAM_OPTIONAL)
        message(FATAL_ERROR "Need ${APP} in a sibling directory.")
      endif()
    endif()
    if(PARAM_DIR AND APP_DIR)
      foreach(DIR ${PARAM_DIR})
        if (EXISTS ${APP_DIR}/${DIR})
          set(MATCH_DIR ${DIR})
        endif()
      endforeach()
      if(NOT MATCH_DIR)
        if(NOT PARAM_OPTIONAL)
          message(FATAL_ERROR "Could not find subdirectory in application.")
        endif()
      else()
        string(APPEND APP_DIR "/${MATCH_DIR}")
      endif()
    elseif(NOT PARAM_DIR)
      get_filename_component(MATCH_DIR ${APP_DIR} NAME)
    endif()
    if(MATCH_DIR)
      message(STATUS "Including ${MATCH_DIR} for project ${CMAKE_PROJECT_NAME}")
      add_subdirectory(${APP_DIR} ${PARAM_TARGET})
      set(${PARAM_TARGET}_DIRECTORY ${APP_DIR} PARENT_SCOPE)
    endif()
  else()
    get_target_property(DIR ${PARAM_TARGET} SOURCE_DIR)
    set(${PARAM_TARGET}_DIRECTORY ${DIR} PARENT_SCOPE)
  endif()
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
endfunction()
