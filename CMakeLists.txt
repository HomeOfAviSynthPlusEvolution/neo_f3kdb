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

set_property(TARGET ${LIBRARY_NAME} PROPERTY CXX_STANDARD 17)
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

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_FLAGS "-march=penryn")
  else()
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_FLAGS "/arch:SSE2")
  endif()
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_FLAGS "/arch:AVX2")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_FLAGS "/arch:AVX512")

  if(CMAKE_GENERATOR_TOOLSET MATCHES "v[0-9]*_xp")
    target_compile_definitions(${LIBRARY_NAME} PRIVATE WINVER=0x502 _WIN32_WINNT=0x502)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:threadSafeInit-")
  endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_FLAGS "/arch:SSE2")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_FLAGS "/arch:AVX2")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_FLAGS "/arch:AVX512")
  target_link_libraries(${LIBRARY_NAME} libmmds)
  target_compile_options(${LIBRARY_NAME} PRIVATE /O3)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_FLAGS "-msse4.1")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_FLAGS "-mavx2;-mfma")
  set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_FLAGS "-mavx512f;-mavx512bw;-mavx512dq;-mavx512vl;-mavx512cd;-mfma")
  target_link_libraries(${LIBRARY_NAME} libmmds)
  target_compile_options(${LIBRARY_NAME} PRIVATE -O3)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    target_compile_options(${LIBRARY_NAME} PRIVATE -O3 -funroll-loops)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
    target_compile_options(${LIBRARY_NAME} PRIVATE -O3 -funroll-loops)
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_FLAGS "-msse4.1")
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_FLAGS "-mavx2;-mfma")
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_FLAGS "-mavx512f;-mavx512bw;-mavx512dq;-mavx512vl;-mavx512cd;-mfma")
  else()
    target_compile_options(${LIBRARY_NAME} PRIVATE -O3 -funroll-loops)
  endif()

else()
  target_compile_options(${LIBRARY_NAME} PRIVATE -O3)
endif()

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(execution HAS_EXECUTION)
if(HAS_EXECUTION)
  add_definitions(-DHAS_EXECUTION)
endif()
if(ENABLE_PAR AND HAS_EXECUTION)
  add_definitions(-DENABLE_PAR)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(${LIBRARY_NAME} tbb)
  endif()
endif()

add_custom_command(
  TARGET ${LIBRARY_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${LIBRARY_NAME}> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:${LIBRARY_NAME}>"
)
