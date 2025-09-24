# Add the appropriate target to generate doxy for used build configuration
#  Adds target to handle doxy, and its installation
# Single-valued parameters:
#   DOX        - Name of .dox file for application
#   TARGET     - Name of application
# Multi-valued parameters:
#   EXTRA_DIRS - Extra directories to add to doxy
#     Used for external directories and optional dependencies
function(ifem_add_doc_target)
  find_package(Doxygen QUIET)
  if(NOT Doxygen_FOUND)
    return()
  endif()

  set(oneValueArgs TARGET DOX)
  set(multiValueArgs EXTRA_DIRS)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT PARAM_DOX)
    set(PARAM_DOX ${PARAM_TARGET})
  endif()

  set(EXTRA_DOXY_PATHS)
  list(APPEND PARAM_EXTRA_DIRS ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
  foreach(DIR ${PARAM_EXTRA_DIRS})
    string(APPEND EXTRA_DOXY_PATHS "${DIR} \\\n")
  endforeach()

  if(IFEM_AS_SUBMODULE)
    configure_file(${IFEM_PATH}/doc/Doxyfile.in Doxyfile)
  else()
    configure_file(doc/Doxyfile.in Doxyfile)
  endif()
  configure_file(doc/${PARAM_DOX}.dox.in ${PARAM_DOX}.dox)

  if(IFEM_INSTALL_DOXY)
    if(NOT CMAKE_INSTALL_DOCDIR)
      set(CMAKE_INSTALL_DOCDIR share/doc)
    endif()
    install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")")
    if(appname STREQUAL "ifem")
      set(dest_path ${CMAKE_INSTALL_DOCDIR}/IFEM)
    else()
      set(dest_path ${CMAKE_INSTALL_DOCDIR}/IFEM/Apps/${PARAM_TARGET})
    endif()
    install(
      DIRECTORY
        ${PROJECT_BINARY_DIR}/doc/html
      DESTINATION
        ${dest_path}
      PATTERN
        *.md5 EXCLUDE
      PATTERN
        *.map EXCLUDE)
  endif(IFEM_INSTALL_DOXY)

  if(IFEM_COMMON_APP_BUILD OR IFEM_AS_SUBMODULE)
    if(NOT TARGET doc)
      add_custom_target(doc)
    endif()
    add_custom_target(
      ${PARAM_TARGET}_doc
      COMMAND
        Doxygen::doxygen ${PROJECT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}
      COMMENT
        "Generating API documentation" VERBATIM
    )
    add_dependencies(doc ${PARAM_APP}_doc)
  else()
    add_custom_target(
      doc
      COMMAND
        Doxygen::doxygen ${PROJECT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}
      COMMENT
        "Generating API documentation" VERBATIM
    )
  endif()
endfunction()
