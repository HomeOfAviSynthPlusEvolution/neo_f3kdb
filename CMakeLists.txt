cmake_minimum_required(VERSION 3.24)

if(POLICY CMP0091)
  cmake_policy(SET CMP0091 NEW)
endif()

project(neo_f3kdb VERSION 10.0.0 LANGUAGES CXX)

option(F3KDB_BUILD_TEST_TOOLS "Build Neo-F3KDB baseline test tools" ON)
option(F3KDB_ENABLE_CLANG_TIDY "Enable clang-tidy static analysis during compilation" OFF)

set(F3KDB_AVISYNTHPLUS_TAG "751ce2021cab1c8405d04825ac2f7e8ac93e8fa0" CACHE STRING "AviSynthPlus tag used when pkg-config cannot find AviSynth headers")
set(F3KDB_VAPOURSYNTH_TAG "R72" CACHE STRING "VapourSynth tag used when pkg-config cannot find VapourSynth headers")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(F3KDB_ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    message(STATUS "Found clang-tidy: ${CLANG_TIDY_EXE}")
    set(F3KDB_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=-Wno-unknown-argument;--extra-arg=-Wno-unused-command-line-argument")
  else()
    message(WARNING "clang-tidy requested but not found in PATH")
  endif()
endif()

if(MSVC AND NOT CMAKE_MSVC_RUNTIME_LIBRARY)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

include(FetchContent)

# 1. Fetch Google Highway
FetchContent_Declare(
  highway
  GIT_REPOSITORY https://github.com/google/highway.git
  GIT_TAG 1.4.0
)
set(HWY_ENABLE_CONTRIB OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(highway)

# 2. Integrate DualSynth2 from Git
set(DS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(DS_BUILD_ACCEPTANCE_PLUGIN OFF CACHE BOOL "" FORCE)
if(MSVC)
  set(DS_MSVC_RUNTIME_LIBRARY "${CMAKE_MSVC_RUNTIME_LIBRARY}" CACHE STRING "" FORCE)
endif()
FetchContent_Declare(
  dualsynth2
  GIT_REPOSITORY https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2.git
  GIT_TAG 0.2
)
FetchContent_MakeAvailable(dualsynth2)

# 3. AviSynth / VapourSynth Headers Configuration
function(f3kdb_fetch_avisynth_headers out_var)
  FetchContent_Declare(
    avisynthplus_headers
    GIT_REPOSITORY https://github.com/AviSynth/AviSynthPlus.git
    GIT_TAG "${F3KDB_AVISYNTHPLUS_TAG}"
    SOURCE_SUBDIR f3kdb-no-cmake
  )
  FetchContent_MakeAvailable(avisynthplus_headers)
  FetchContent_GetProperties(avisynthplus_headers)
  set(${out_var} "${avisynthplus_headers_SOURCE_DIR}/avs_core/include" PARENT_SCOPE)
endfunction()

function(f3kdb_fetch_vapoursynth_headers out_var)
  FetchContent_Declare(
    vapoursynth_headers
    GIT_REPOSITORY https://github.com/vapoursynth/vapoursynth.git
    GIT_TAG "${F3KDB_VAPOURSYNTH_TAG}"
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR f3kdb-no-cmake
  )
  FetchContent_MakeAvailable(vapoursynth_headers)
  FetchContent_GetProperties(vapoursynth_headers)
  set(${out_var} "${vapoursynth_headers_SOURCE_DIR}/include" PARENT_SCOPE)
endfunction()

function(f3kdb_add_vapoursynth_compat_headers target)
  set(compat_dir "${CMAKE_CURRENT_BINARY_DIR}/generated/vapoursynth-compat")
  file(MAKE_DIRECTORY "${compat_dir}/vapoursynth")
  foreach(header IN ITEMS
    VSConstants4.h
    VSHelper.h
    VSHelper4.h
    VSScript.h
    VSScript4.h
    VapourSynth.h
    VapourSynth4.h
  )
    file(WRITE
      "${compat_dir}/vapoursynth/${header}"
      "#pragma once\n#include <${header}>\n"
    )
  endforeach()
  target_include_directories(${target} INTERFACE "${compat_dir}")
endfunction()

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(AVISYNTH QUIET avisynth)
  pkg_check_modules(VAPOURSYNTH QUIET vapoursynth)
endif()

add_library(f3kdb_avisynth_headers INTERFACE)
if(AVISYNTH_FOUND)
  target_include_directories(f3kdb_avisynth_headers INTERFACE ${AVISYNTH_INCLUDE_DIRS})
  target_compile_options(f3kdb_avisynth_headers INTERFACE ${AVISYNTH_CFLAGS_OTHER})
  message(STATUS "Using AviSynth headers from pkg-config")
else()
  f3kdb_fetch_avisynth_headers(F3KDB_AVISYNTH_INCLUDE_DIRS)
  target_include_directories(f3kdb_avisynth_headers INTERFACE ${F3KDB_AVISYNTH_INCLUDE_DIRS})
  message(STATUS "Using AviSynth headers from FetchContent (${F3KDB_AVISYNTHPLUS_TAG})")
endif()

add_library(f3kdb_vapoursynth_headers INTERFACE)
if(VAPOURSYNTH_FOUND)
  target_include_directories(f3kdb_vapoursynth_headers INTERFACE ${VAPOURSYNTH_INCLUDE_DIRS})
  target_compile_options(f3kdb_vapoursynth_headers INTERFACE ${VAPOURSYNTH_CFLAGS_OTHER})
  message(STATUS "Using VapourSynth headers from pkg-config")
else()
  f3kdb_fetch_vapoursynth_headers(F3KDB_VAPOURSYNTH_INCLUDE_DIRS)
  target_include_directories(f3kdb_vapoursynth_headers INTERFACE ${F3KDB_VAPOURSYNTH_INCLUDE_DIRS})
  message(STATUS "Using VapourSynth headers from FetchContent (${F3KDB_VAPOURSYNTH_TAG})")
endif()
f3kdb_add_vapoursynth_compat_headers(f3kdb_vapoursynth_headers)

# 4. Version and Resource Generation
string(REGEX MATCH "^([0-9.]+)" CORE_VERSION_STRING "${PROJECT_VERSION}")
if (CORE_VERSION_STRING)
    if(PROJECT_VERSION MATCHES "-")
        set(IS_PRERELEASE 1)
    else()
        set(IS_PRERELEASE 0)
    endif()
else()
    set(CORE_VERSION_STRING "10.0.0")
    set(IS_PRERELEASE 0)
endif()

string(REPLACE "." ";" VERSION_NUMERIC_PARTS_LIST "${CORE_VERSION_STRING}")
list(LENGTH VERSION_NUMERIC_PARTS_LIST NUM_PARTS)

set(VERSION_MAJOR 0)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION_BUILD 0)

if(NUM_PARTS GREATER_EQUAL 1)
    list(GET VERSION_NUMERIC_PARTS_LIST 0 VERSION_MAJOR)
endif()
if(NUM_PARTS GREATER_EQUAL 2)
    list(GET VERSION_NUMERIC_PARTS_LIST 1 VERSION_MINOR)
endif()
if(NUM_PARTS GREATER_EQUAL 3)
    list(GET VERSION_NUMERIC_PARTS_LIST 2 VERSION_PATCH)
endif()
if(NUM_PARTS GREATER_EQUAL 4)
    list(GET VERSION_NUMERIC_PARTS_LIST 3 VERSION_BUILD)
endif()

set(PROJECT_VERSION_STRING_FULL "v${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

set(F3KDB_GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
configure_file(
  "${PROJECT_SOURCE_DIR}/src/version.hpp.in"
  "${F3KDB_GENERATED_DIR}/version.hpp"
)

if (WIN32)
    set(FILE_DESCRIPTION "Neo-F3KDB")
    set(INTERNAL_NAME "Neo-F3KDB")
    set(ORIGINAL_FILENAME "neo-f3kdb.dll")
    set(PRODUCT_NAME "Neo-F3KDB")

    set(RC_FILEFLAGS_LINE "FILEFLAGS   0x0L")
    set(RC_FILE_OUT "${F3KDB_GENERATED_DIR}/version.rc")
    configure_file(
        "${PROJECT_SOURCE_DIR}/src/version.rc.in"
        "${RC_FILE_OUT}"
        @ONLY
    )
endif()

# 5. Target Sources Setup
set(LIBRARY_NAME "neo-f3kdb")

set(F3KDB_SOURCES
  src/plugin/plugin_entry.cpp
  src/plugin/f3kdb_filter.cpp
  src/plugin/f3kdb_filter.hpp
  src/bit_utils.h
  src/compiler_compat.h
  src/constants.h
  src/core.cpp
  src/core.h
  src/cpuid.cpp
  src/dither_high.h
  src/f3kdb.h
  src/f3kdb.hpp
  src/flash3kyuu_deband_impl_c.cpp
  src/impl_dispatch.cpp
  src/impl_dispatch.h
  src/impl_dispatch_decl.h
  src/neo_f3kdb.hpp
  src/pixel_proc_c.h
  src/pixel_proc_c_16bit.h
  src/pixel_proc_c_high_bit_depth_common.h
  src/pixel_proc_c_high_f_s_dithering.h
  src/pixel_proc_c_high_no_dithering.h
  src/pixel_proc_c_high_ordered_dithering.h
  src/process_plane_context.cpp
  src/process_plane_context.h
  src/random.cpp
  src/random.h
  src/sse2neon.h
  src/sse_utils.h
  src/utils.h
)

if(WIN32)
  list(APPEND F3KDB_SOURCES "${RC_FILE_OUT}")
endif()

add_library(${LIBRARY_NAME} SHARED ${F3KDB_SOURCES})
set_property(TARGET ${LIBRARY_NAME} PROPERTY CXX_STANDARD 23)

target_include_directories(${LIBRARY_NAME}
  PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/include/dualsynth"
    "${F3KDB_GENERATED_DIR}"
)

# GNU-style optimization flags
if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
    target_compile_options(${LIBRARY_NAME} PRIVATE -O3 -funroll-loops)
endif()

# Only apply SIMD flags if we are on a capable architecture (x86/x86_64).
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i.86")
    message(STATUS "x86/x86_64 architecture detected. Configuring SIMD instruction sets.")

    target_sources(${LIBRARY_NAME} PRIVATE
        "src/flash3kyuu_deband_avx2_base.h"
        "src/flash3kyuu_deband_avx512_base.h"
        "src/flash3kyuu_deband_impl_avx2.cpp"
        "src/flash3kyuu_deband_impl_avx512.cpp"
        "src/flash3kyuu_deband_impl_sse4.cpp"
        "src/flash3kyuu_deband_sse_base.h"
    )

    target_include_directories(${LIBRARY_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/VCL2")

    set_source_files_properties("src/flash3kyuu_deband_impl_sse4.cpp" PROPERTIES COMPILE_OPTIONS
        "$<$<CXX_COMPILER_ID:MSVC>:/arch:SSE2>$<$<AND:$<CXX_COMPILER_ID:Intel>,$<PLATFORM_ID:Windows>>:/arch:SSE4.1>$<$<NOT:$<OR:$<CXX_COMPILER_ID:MSVC>,$<AND:$<CXX_COMPILER_ID:Intel>,$<PLATFORM_ID:Windows>>>>:-msse4.1>"
    )

    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        set_source_files_properties("src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
    elseif(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
        set_source_files_properties("src/flash3kyuu_deband_impl_avx2.cpp" PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
    endif()

    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        set_source_files_properties("src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_OPTIONS "/arch:AVX512")
    elseif(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
        set_source_files_properties("src/flash3kyuu_deband_impl_avx512.cpp" PROPERTIES COMPILE_OPTIONS "-mavx512f;-mavx512bw;-mavx512dq;-mavx512vl;-mavx512cd;-mfma")
    endif()

    target_compile_definitions(${LIBRARY_NAME} PRIVATE HAS_SSE4 HAS_AVX2 HAS_AVX512)
else()
    message(STATUS "Non-x86 architecture detected (${CMAKE_SYSTEM_PROCESSOR}). Skipping some SIMD-specific source files.")

    target_include_directories(${LIBRARY_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/src/arm")
    target_sources(${LIBRARY_NAME} PRIVATE "src/flash3kyuu_deband_impl_sse4.cpp")
    target_compile_definitions(${LIBRARY_NAME} PRIVATE HAS_SSE4 INSTRSET=6)

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        target_compile_options(${LIBRARY_NAME} PRIVATE -march=armv8-a+simd -fsigned-char)
    else()
        target_compile_options(${LIBRARY_NAME} PRIVATE -mfpu=neon -fsigned-char)
    endif()
endif()

target_link_libraries(${LIBRARY_NAME} PRIVATE "$<$<OR:$<CXX_COMPILER_ID:Intel>,$<CXX_COMPILER_ID:IntelLLVM>>:libmmds>")

# Link Google Highway, Custom Headers and DualSynth2
target_link_libraries(${LIBRARY_NAME} PRIVATE hwy f3kdb_avisynth_headers f3kdb_vapoursynth_headers)
target_link_libraries(${LIBRARY_NAME} PRIVATE DualSynth::dualsynth)

if(CMAKE_DL_LIBS)
  target_link_libraries(${LIBRARY_NAME} PRIVATE ${CMAKE_DL_LIBS})
endif()

if(F3KDB_CLANG_TIDY)
  set_target_properties(${LIBRARY_NAME} PROPERTIES CXX_CLANG_TIDY "${F3KDB_CLANG_TIDY}")
endif()

# 6. Build Baseline Test Tools
if(F3KDB_BUILD_TEST_TOOLS)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  add_executable(f3kdb_avs_frame_dump
    tests/tools/avs_frame_dump.cpp
  )
  set_property(TARGET f3kdb_avs_frame_dump PROPERTY CXX_STANDARD 17)
  target_link_libraries(f3kdb_avs_frame_dump PRIVATE f3kdb_avisynth_headers)
  if(NOT WIN32)
    target_link_libraries(f3kdb_avs_frame_dump PRIVATE dl)
  endif()

  enable_testing()

  # Register the baseline verification tasks with CTest
  foreach(backend IN ITEMS purec)
    foreach(tier IN ITEMS smoke compat)
      add_test(
        NAME baseline_${backend}_${tier}
        COMMAND
          ${Python3_EXECUTABLE} -m tests.baseline.baseline_runner verify
          --plugin $<TARGET_FILE:${LIBRARY_NAME}>
          --avs-dump $<TARGET_FILE:f3kdb_avs_frame_dump>
          --golden "${PROJECT_SOURCE_DIR}/tests/baseline/golden/${tier}.json"
          --tier ${tier}
          --backend ${backend}
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
      )
    endforeach()
  endforeach()
endif()
