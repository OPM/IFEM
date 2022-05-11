IF(VTFWRITER_INCLUDES AND VTFWRITER_LIBRARIES)
  SET(VTFWRITER_FIND_QUIETLY TRUE)
ENDIF(VTFWRITER_INCLUDES AND VTFWRITER_LIBRARIES)

UNSET(VTFWRITER_INCLUDES CACHE)
UNSET(VTFWRITER_LIBRARIES CACHE)

IF (NOT VTFAPI OR VTFAPI GREATER 1)
  FIND_PATH(VTFWRITER_INCLUDES
    NAMES VTFXAPI.h
    PATHS $ENV{HOME}/include
    /usr/local/include
    /sima/libs/VTFx/include
  )
  IF(VTFWRITER_INCLUDES)
    SET(VTFAPI 2)
  ENDIF(VTFWRITER_INCLUDES)
ENDIF (NOT VTFAPI OR VTFAPI GREATER 1)

IF (NOT VTFWRITER_INCLUDES)
  FIND_PATH(VTFWRITER_INCLUDES
    NAMES VTFAPI.h
    PATHS $ENV{HOME}/include
    /usr/local/include
    /sima/libs/GLviewExpressWriter/include
  )
  IF(VTFWRITER_INCLUDES)
    SET(VTFAPI 1)
  ENDIF(VTFWRITER_INCLUDES)
ENDIF(NOT VTFWRITER_INCLUDES)

IF(VTFWRITER_INCLUDES)
  IF (${VTFAPI} LESS 2)
    FIND_LIBRARY(VTFWRITER_LIBRARIES
      NAMES VTFExpressAPI
      PATHS $ENV{HOME}/lib
      /usr/local/lib
      /sima/libs/GLviewExpressWriter/lib
    )
    MESSAGE(STATUS "Using VTF libraries")
  ELSE(${VTFAPI} LESS 2)
    FIND_LIBRARY(VTFWRITER_LIBRARIES
      NAMES GLviewExpressWriter2d
      PATHS $ENV{HOME}/lib
      /usr/local/lib
      /sima/libs/VTFx/lib
    )
    FIND_LIBRARY(ZIP_LIBRARIES
      NAMES ziparch
      PATHS $ENV{HOME}/lib
      /usr/local/lib
      /sima/libs/VTFx/lib
    )
    SET(VTFWRITER_LIBRARIES ${VTFWRITER_LIBRARIES} ${ZIP_LIBRARIES})
    MESSAGE(STATUS "Using VTFx libraries")
  ENDIF(${VTFAPI} LESS 2)
ENDIF(VTFWRITER_INCLUDES)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VTFWriter DEFAULT_MSG
                                  VTFWRITER_INCLUDES VTFWRITER_LIBRARIES)

MARK_AS_ADVANCED(VTFAPI VTFWRITER_INCLUDES VTFWRITER_LIBRARIES)
