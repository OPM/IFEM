# - Find TinyXML2
# Find the native TinyXML includes and library
#
#   TinyXML2_FOUND       - True if TinyXML found.
# Adds target tinyxml2::tinyxml2

find_path(TinyXML2_INCLUDE_DIR "tinyxml2.h"
          PATH_SUFFIXES "tinyxml2" )

find_library(TinyXML2_LIBRARIES
             NAMES "tinyxml2"
             PATH_SUFFIXES "tinyxml2" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( "TinyXML2" DEFAULT_MSG TinyXML2_INCLUDE_DIR TinyXML2_LIBRARIES )

add_library(tinyxml2::tinyxml2 IMPORTED UNKNOWN)
set_target_properties(tinyxml2::tinyxml2 PROPERTIES
                      IMPORTED_LOCATION ${TinyXML2_LIBRARIES}
                      INCLUDE_DIRECTORIES ${TinyXML2_INCLUDE_DIR})

mark_as_advanced(TinyXML2_INCLUDE_DIR TinyXML_LIBRARIES)
