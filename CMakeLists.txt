cmake_minimum_required(VERSION 3.15)
project(neo_f3kdb CXX)

set(LIBRARY_NAME "neo-f3kdb")
add_library(${LIBRARY_NAME} SHARED)

target_sources(${LIBRARY_NAME} PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/bit_utils.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/compiler_compat.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/constants.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cpuid.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/dither_high.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/f3kdb.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/f3kdb.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_avx2_base.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_avx512_base.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_c.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_sse_base.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/impl_dispatch.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/impl_dispatch.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/impl_dispatch_decl.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/neo_f3kdb.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pixel_proc_c.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pixel_proc_c_16bit.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pixel_proc_c_high_bit_depth_common.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pixel_proc_c_high_f_s_dithering.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pixel_proc_c_high_no_dithering.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pixel_proc_c_high_ordered_dithering.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/process_plane_context.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/process_plane_context.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/random.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/random.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/sse2neon.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/sse_utils.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/utils.h"
)

target_compile_features(${LIBRARY_NAME} PRIVATE cxx_std_17)
option(ENABLE_PAR "Enable C++17 Parallel Execution" ON)

find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} describe --first-parent --tags --always OUTPUT_VARIABLE GIT_REPO_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX REPLACE "(r[0-9]+).*" "\\1" VERSION ${GIT_REPO_VERSION})

configure_file(
  "${PROJECT_SOURCE_DIR}/src/version.hpp.in"
  "${PROJECT_SOURCE_DIR}/src/version.hpp"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/version.rc.in"
  "${PROJECT_SOURCE_DIR}/src/version.rc"
)

if (NOT MSVC)
  find_package(PkgConfig REQUIRED)

  pkg_check_modules(AVISYNTH avisynth)
  if(AVISYNTH_FOUND)
    target_include_directories(${LIBRARY_NAME} PRIVATE ${AVISYNTH_INCLUDE_DIRS})
  else()
    target_include_directories(${LIBRARY_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/avisynth")
  endif()

  pkg_check_modules(VAPOURSYNTH vapoursynth)
  if(VAPOURSYNTH_FOUND)
    target_include_directories(${LIBRARY_NAME} PRIVATE ${VAPOURSYNTH_INCLUDE_DIRS})
  else()
    target_include_directories(${LIBRARY_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/vapoursynth")
  endif()
else()
  target_include_directories(${LIBRARY_NAME} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include/avisynth"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/vapoursynth"
  )
endif()

target_include_directories(${LIBRARY_NAME} PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/dualsynth"
  "${CMAKE_CURRENT_SOURCE_DIR}/VCL2"
)

# We apply aggressive flags ONLY to compilers with a GNU-style frontend.
if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
    target_compile_options(${LIBRARY_NAME} PRIVATE -O3 -funroll-loops)
endif()

# Only apply SIMD flags if we are on a capable architecture (x86/x86_64).
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i.86")
    message(STATUS "x86/x86_64 architecture detected. Configuring SIMD instruction sets.")

    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_OPTIONS
        "$<$<CXX_COMPILER_ID:MSVC>:/arch:SSE2>$<$<AND:$<CXX_COMPILER_ID:Intel>,$<PLATFORM_ID:Windows>>:/arch:SSE4.1>$<$<NOT:$<OR:$<CXX_COMPILER_ID:MSVC>,$<AND:$<CXX_COMPILER_ID:Intel>,$<PLATFORM_ID:Windows>>>>:-msse4.1>"
    )

    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
    elseif(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
        set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
    endif()

    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_OPTIONS "/arch:AVX512")
    elseif(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
        set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_OPTIONS "-mavx512f;-mavx512bw;-mavx512dq;-mavx512vl;-mavx512cd;-mfma")
    endif()
else()
    message(STATUS "Non-x86 architecture detected (${CMAKE_SYSTEM_PROCESSOR}). Skipping SIMD-specific source files.")
endif()

target_link_libraries(${LIBRARY_NAME} PRIVATE "$<$<OR:$<CXX_COMPILER_ID:Intel>,$<CXX_COMPILER_ID:IntelLLVM>>:libmmds>")

# Handle legacy Windows XP support if the specific toolset is used
if(CMAKE_GENERATOR_TOOLSET MATCHES "v[0-9]+_xp")
    message(STATUS "Windows XP toolset detected. Applying compatibility settings.")
    target_compile_definitions(${LIBRARY_NAME} PRIVATE WINVER=0x502 _WIN32_WINNT=0x502)
    target_compile_options(${LIBRARY_NAME} PRIVATE "$<$<CXX_COMPILER_ID:MSVC>:/Zc:threadSafeInit->")
endif()

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(execution HAS_EXECUTION)
if(HAS_EXECUTION)
  target_compile_definitions(${LIBRARY_NAME} PRIVATE HAS_EXECUTION)
endif()
if(ENABLE_PAR AND HAS_EXECUTION)
  target_compile_definitions(${LIBRARY_NAME} PRIVATE ENABLE_PAR)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(${LIBRARY_NAME} PRIVATE tbb)
  endif()
endif()

add_custom_command(
  TARGET ${LIBRARY_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${LIBRARY_NAME}> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:${LIBRARY_NAME}>"
)
