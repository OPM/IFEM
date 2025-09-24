find_library(
  SPR8_LIBRARIES
  NAMES
   SPR_I8
  PATHS
    $ENV{HOME}/lib
    $ENV{HOME}/.local/lib
)

find_library(
  SPR4_LIBRARIES
  NAMES
    SPR
  PATHS
    $ENV{HOME}/lib
    $ENV{HOME}/.local/lib
)

mark_as_advanced(SPR8_LIBRARIES SPR4_LIBRARIES)

if (SPR4_LIBRARIES OR SPR8_LIBRARIES)
  find_library(
    ASM_LIBRARIES
    NAMES
      ASM SAM
    PATHS
      $ENV{HOME}/lib
      $ENV{HOME}/.local/lib
  )
  mark_as_advanced(ASM_LIBRARIES)
  if(ASM_LIBRARIES)
    add_library(SPR::ASM UNKNOWN IMPORTED)
    set_target_properties(SPR::ASM PROPERTIES IMPORTED_LOCATION ${ASM_LIBRARIES})
    target_compile_definitions(SPR::ASM INTERFACE USE_F77SAM=1)
  endif()
endif()

if(SPR4_LIBRARIES AND ASM_LIBRARIES)
  add_library(SPR::INT4 UNKNOWN IMPORTED)
  set_target_properties(SPR::INT4 PROPERTIES IMPORTED_LOCATION ${SPR4_LIBRARIES})
  if(TARGET SPR::ASM)
    target_link_libraries(SPR::INT4 INTERFACE SPR::ASM)
  endif()
endif()

if(SPR8_LIBRARIES)
  add_library(SPR::INT8 UNKNOWN IMPORTED)
  set_target_properties(SPR::INT8 PROPERTIES IMPORTED_LOCATION ${SPR8_LIBRARIES})
  if(TARGET SPR::ASM)
    target_link_libraries(SPR::INT8 INTERFACE SPR::ASM)
  endif()
  target_compile_definitions(SPR::INT8 INTERFACE USE_INT64=1)
endif()

if (SPR_USE_INT8)
  set(SPR_LIBRARIES ${SPR4_LIBRARIES})
else()
  set(SPR_LIBRARIES ${SPR4_LIBRARIES})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPR DEFAULT_MSG SPR_LIBRARIES)

if(SPR_FOUND)
  add_library(SPR::SPR UNKNOWN IMPORTED)
  set_target_properties(SPR::SPR PROPERTIES IMPORTED_LOCATION ${SPR_LIBRARIES})
  if (SPR_USE_INT8)
    target_link_libraries(SPR::SPR INTERFACE SPR::INT8)
  else()
    target_link_libraries(SPR::SPR INTERFACE SPR::INT4)
  endif()
  target_compile_definitions(SPR::SPR INTERFACE HAS_SPR=1)
endif()
