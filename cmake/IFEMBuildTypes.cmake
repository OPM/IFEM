# Set default build type if none is set and handle some build type aliases

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Custom profiles
if(${CMAKE_BUILD_TYPE} MATCHES "Nopt")
  set(CMAKE_BUILD_TYPE Debug)
elseif(${CMAKE_BUILD_TYPE} MATCHES "Nomp")
  set(CMAKE_BUILD_TYPE Release)
endif()
