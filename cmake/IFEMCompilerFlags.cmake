# Probe and add some compiler flags for warnings and optimization

include(CheckCXXCompilerFlag)

if(CMAKE_BUILD_TYPE MATCHES "Release" AND IFEM_USE_NATIVE)
  check_cxx_compiler_flag("-mtune=native" HAVE_MTUNE)
  if(HAVE_MTUNE)
    target_compile_options(IFEM PUBLIC -mtune=native)
  endif(HAVE_MTUNE)
endif()

check_cxx_compiler_flag(-Wall HAS_WALL)
if(HAS_WALL)
  target_compile_options(IFEM PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
endif()

check_cxx_compiler_flag(-Wno-parentheses HAS_PARENTHESES)
if(HAS_PARENTHESES)
  target_compile_options(IFEM PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-Wno-parentheses>)
endif()
