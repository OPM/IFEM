include(CheckIPOSupported)
include(TestCXXAcceptsFlag)

# ------------------------------------------------------------------------------
# ifem_lto
#
# Configure link-time optimization (IPO) for a target.
#
# This function enables LTO (Link Time Optimization) on the specified target.
# It has extended interfaces for GCC and Clang settings.
#
# Enablement is controlled using the IFEM_USE_LTO option.
#
# Incremental cache directories are automatically managed under:
#   ${CMAKE_BINARY_DIR}/LTOCache/<config>
#
# Usage:
#
#   ifem_ipo(
#     TARGET <target>
#     [CONFIGURATION_TYPES <config1;config2;...>]
#   )
#
# Arguments:
#   TARGET               (required) Target to which optimization will be applied.
#   CONFIGURATION_TYPES  (optional) List of configuration names (e.g. Release).
#                        If omitted, the default list from:
#                          IFEM_LTO_CONFIGURATION_TYPES
#                        is used.
#
# Cache Variables:
#   IFEM_LTO_JOBS                 (e.g. 1, 8)
#   IFEM_LTO_CACHE_POLICY         (e.g. "prune_after=604800")
#   IFEM_LTO_CONFIGURATION_TYPES  (e.g. "Release;RelWithDebInfo")
#
# ------------------------------------------------------------------------------

set(IFEM_LTO_CONFIGURATION_TYPES "Release;RelWithDebInfo" CACHE
    STRING "Default configuration types to apply interprocedural optimization.")
set(IFEM_LTO_JOBS 1 CACHE STRING "Default interprocedural optimization jobs. Set to 'auto' to use all cores")
set(IFEM_LTO_CACHE_POLICY "" CACHE STRING "Default interprocedural optimization cache policy.")
mark_as_advanced(IFEM_LTO_CACHE_POLICY)
mark_as_advanced(IFEM_LTO_CONFIGURATION_TYPES)

# define possible values for LTO (order matters: first one in the list is used as default type)
set(ifem_ipo_types)
check_ipo_supported(LANGUAGES C CXX RESULT ipo_supported)

if(IFEM_USE_LTO AND NOT ipo_supported)
  message(STATUS "Link-time optimization not supported")
  set(IFEM_USE_LTO OFF)
endif()

if(IFEM_USE_LTO)
  message(STATUS "Link-time optimization enabled")
endif()

function(ifem_lto)
  if(NOT IFEM_USE_LTO)
    return()
  endif()

  cmake_parse_arguments(PARAM "" "TARGET" "CONFIGURATION_TYPES" ${ARGN})

  # if not set, use the default configuration types
  if(NOT PARAM_CONFIGURATION_TYPES)
    set(PARAM_CONFIGURATION_TYPES ${IFEM_LTO_CONFIGURATION_TYPES})
  endif()

  foreach(config ${PARAM_CONFIGURATION_TYPES})
    if(CMAKE_CONFIGURATION_TYPES MATCHES ${config} OR CMAKE_BUILD_TYPE MATCHES ${config})
      set(LTO_CACHE_PATH ${CMAKE_BINARY_DIR}/LTOCache/${config})

      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gcc/Optimize-Options.html
        target_compile_options(${PARAM_TARGET} PRIVATE $<$<CONFIG:${config}>:-flto=${IFEM_LTO_JOBS}>)
        target_link_options(${PARAM_TARGET} INTERFACE $<$<CONFIG:${config}>:-flto=${IFEM_LTO_JOBS}>)

        if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 15.0.0)
          # use incremental LTO https://gcc.gnu.org/onlinedocs/gcc-15.2.0/gcc/Optimize-Options.html
          target_compile_options(${PARAM_TARGET} PRIVATE $<$<CONFIG:${config}>:-flto-incremental=${LTO_CACHE_PATH}>)
          target_link_options(${PARAM_TARGET} INTERFACE $<$<CONFIG:${config}>:-flto-incremental=${LTO_CACHE_PATH}>)

          # Configure cache directory
          file(MAKE_DIRECTORY ${LTO_CACHE_PATH})
          set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${LTO_CACHE_PATH})
        endif()
      elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        # Use ThinLTO and reuse cache (see https://releases.llvm.org/20.1.0/tools/clang/docs/ThinLTO.html)

        if(${CMAKE_CXX_COMPILER_LINKER_ID} MATCHES "GNU")
          message(FATAL_ERROR "Cannot use ThinLTO with BFD linker")
        endif()

        target_compile_options(${PARAM_TARGET} PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:${config}>>:-flto=thin>)
        target_link_options(${PARAM_TARGET} INTERFACE $<$<CONFIG:${config}>:-flto=thin>)

        # Configure cache directory
        file(MAKE_DIRECTORY ${LTO_CACHE_PATH})
        set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${LTO_CACHE_PATH})

        # Parallel and cache options are linker dependent
        if(${CMAKE_CXX_COMPILER_LINKER_ID} MATCHES "LLD")
          target_link_options(${PARAM_TARGET} INTERFACE
            $<$<CONFIG:${config}>:LINKER:--thinlto-jobs=${IFEM_LTO_JOBS}>
            $<$<CONFIG:${config}>:LINKER:--thinlto-cache-dir=${LTO_CACHE_PATH}>
          )
          if(IFEM_LTO_CACHE_POLICY)
            target_link_options(${PARAM_TARGET} INTERFACE
              $<$<CONFIG:${config}>:LINKER:--thinlto-cache-policy=${IFEM_LTO_CACHE_POLICY}>
            )
          endif()
        elseif(${CMAKE_CXX_COMPILER_LINKER_ID} MATCHES "GNU|GNUgold|MOLD")
          target_link_options(${PARAM_TARGET} INTERFACE
            $<$<CONFIG:${config}>:LINKER:-plugin-opt,jobs=${IFEM_LTO_JOBS}>
            $<$<CONFIG:${config}>:LINKER:-plugin-opt,cache-dir=${LTO_CACHE_PATH}>
          )
          if(IFEM_LTO_CACHE_POLICY)
            target_link_options(${PARAM_TARGET} INTERFACE
              $<$<CONFIG:${config}>:LINKER:-plugin-opt,cache-policy=${IFEM_LTO_CACHE_POLICY}>
            )
          endif()
        elseif(${CMAKE_CXX_COMPILER_LINKER_ID} MATCHES "AppleClang")
          target_link_options(${PARAM_TARGET} INTERFACE
            $<$<CONFIG:${config}>:LINKER:-mllvm,-threads=${IFEM_LTO_JOBS}>
            $<$<CONFIG:${config}>:LINKER:-cache_path_lto,${LTO_CACHE_PATH}>
          )
        else()
          message(DEBUG "LTO for Clang is supported, but linker type '${CMAKE_CXX_COMPILER_LINKER_ID}' is unrecognized. Skip adding parallelism or cache flags.")
        endif()
      endif()
    endif()
  endforeach()
endfunction()
