find_path(
  VTFXWRITER_INCLUDE_DIRS
  NAMES
    VTFXAPI.h
  PATHS
    $ENV{HOME}/include
    $ENV{HOME}/.local/include
)

find_path(
  VTFWRITER_INCLUDE_DIRS
  NAMES
   VTFAPI.h
  PATHS
    $ENV{HOME}/include
    $ENV{HOME}/.local/include
)


find_library(
  VTFWRITER_LIBRARIES
  NAMES
    VTFExpressAPI
  PATHS
    $ENV{HOME}/lib
    $ENV{HOME}/.local/lib
)

find_library(
  VTFXWRITER_LIBRARIES
  NAMES
    GLviewExpressWriter2d
  PATHS
    $ENV{HOME}/lib
    $ENV{HOME}/.local/lib
)

find_library(
  ZIP_LIBRARIES
  NAMES
    ziparch
  PATHS
    $ENV{HOME}/lib
    $ENV{HOME}/.local/lib
)

if(VTFWRITER_INCLUDE_DIRS AND VTFWRITER_LIBRARIES)
  add_library(VTFWriter::VTF UNKNOWN IMPORTED)
  set_target_properties(VTFWriter::VTF PROPERTIES IMPORTED_LOCATION ${VTFWRITER_LIBRARIES})
  target_include_directories(VTFWriter::VTF INTERFACE ${VTFWRITER_INCLUDE_DIRS})
  target_compile_definitions(VTFWriter::VTF INTERFACE HAS_VTFAPI=1)
endif()

if(VTFXWRITER_INCLUDE_DIRS AND VTFXWRITER_LIBRARIES AND ZIP_LIBRARIES)
  add_library(VTFWriter::VTFx UNKNOWN IMPORTED)
  set_target_properties(VTFWriter::VTFx PROPERTIES IMPORTED_LOCATION ${VTFXWRITER_LIBRARIES})
  target_include_directories(VTFWriter::VTFx INTERFACE ${VTFXWRITER_INCLUDE_DIRS})
  target_link_libraries(VTFWriter::VTFx INTERFACE ${ZIP_LIBRARIES})
  target_compile_definitions(VTFWriter::VTFx INTERFACE HAS_VTFAPI=2)
  set(VTFWRITER_INCLUDE_DIRS ${VTFXWRITER_INCLUDE_DIRS})
  set(VTFWRITER_LIBRARIES ${VTFXWRITER_LIBRARIES})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VTFWriter DEFAULT_MSG
                                  VTFWRITER_INCLUDE_DIRS VTFWRITER_LIBRARIES)

mark_as_advanced(VTFWRITER_INCLUDE_DIRS VTFWRITER_LIBRARIES
                 VTFXWRITER_INCLUDE_DIRS VTFXWRITER_LIBRARIES ZIP_LIBRARIES)
