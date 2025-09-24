# Probes for portability issues

include(CheckFunctionExists)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_FLAGS)
check_function_exists(get_current_dir_name unistd.h HAVE_GET_CURRENT_DIR_NAME) # lacks on osx
if(HAVE_GET_CURRENT_DIR_NAME)
  target_compile_definitions(IFEM PUBLIC HAVE_GET_CURRENT_DIR_NAME=1)
endif()
