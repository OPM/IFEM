# Add the appropriate target to generate doxy for used build configuration
macro(add_doc_target appname dox)
  if(IFEM_AS_SUBMODULE OR IFEM_INTREE_BUILD)
    configure_file(${IFEM_PATH}/doc/Doxyfile.in Doxyfile)
  else()
    configure_file(doc/Doxyfile.in Doxyfile)
  endif()
  configure_file(doc/${dox}.dox.in ${dox}.dox)

  if(IFEM_INSTALL_DOXY)
    install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
    install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}/Apps/${appname}
            COMPONENT doc
            PATTERN *.md5 EXCLUDE
            PATTERN *.map EXCLUDE)
  endif(IFEM_INSTALL_DOXY)

  if(IFEM_COMMON_APP_BUILD)
    if(NOT TARGET doc)
      add_custom_target(doc)
    endif()
    add_custom_target(${appname}_doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Generating API documentation" VERBATIM)
    add_dependencies(doc ${appname}_doc)
  else()
    add_custom_target(doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Generating API documentation" VERBATIM)
  endif()
endmacro()
