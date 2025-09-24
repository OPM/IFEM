find_path(
  SuperLU_INCLUDES
  NAMES
    slu_ddefs.h
  PATHS
    $ENV{HOME}/include
  PATH_SUFFIXES
    superlu
)

find_library(
  SuperLU_LIBRARIES
  NAMES
    superlu
  PATHS
    $ENV{HOME}/lib
)

find_path(
  SuperLU_MT_INCLUDES
  NAMES
    slu_mt_ddefs.h
  PATHS
    $ENV{HOME}/include
  PATH_SUFFIXES
    superlu_mt superlu
)
find_library(
  SuperLU_MT_LIBRARIES
  NAMES
    superlu_mt
  PATHS
    $ENV{HOME}/lib
)

include(FindPackageHandleStandardArgs)
if(SuperLU_LIBRARIES)
  find_package_handle_standard_args(SuperLU DEFAULT_MSG
                                    SuperLU_INCLUDES SuperLU_LIBRARIES)
  add_library(SuperLU::SuperLU IMPORTED UNKNOWN)
  set_target_properties(SuperLU::SuperLU PROPERTIES IMPORTED_LOCATION ${SuperLU_LIBRARIES})
  target_include_directories(SuperLU::SuperLU INTERFACE ${SuperLU_INCLUDES})
  target_compile_definitions(SuperLU::SuperLU INTERFACE HAS_SUPERLU=1)
endif()

if(SuperLU_MT_LIBRARIES)
  find_package(Threads REQUIRED)
  find_package_handle_standard_args(SuperLU_MT DEFAULT_MSG
                                    SuperLU_MT_INCLUDES SuperLU_MT_LIBRARIES)
  add_library(SuperLU::SuperLUMT IMPORTED UNKNOWN)
  set_target_properties(SuperLU::SuperLUMT PROPERTIES IMPORTED_LOCATION ${SuperLU_MT_LIBRARIES})
  target_include_directories(SuperLU::SuperLUMT INTERFACE ${SuperLU_MT_INCLUDES})
  target_link_libraries(SuperLU::SuperLUMT INTERFACE ${CMAKE_THREAD_LIBS_INIT})
  target_compile_definitions(SuperLU::SuperLUMT INTERFACE HAS_SUPERLU=1 HAS_SUPERLU_MT=1)
endif()

mark_as_advanced(SuperLU_INCLUDES SuperLU_LIBRARIES SUPERLU_DEFINITIONS
                 SuperLU_MT_INCLUDES SuperLU_MT_LIBRARIES)
