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
  src/constants.h
  src/compiler_compat.h
  src/core/aligned_buffer.hpp
  src/core/dither_floyd_steinberg.hpp
  src/core/dither_info.cpp
  src/core/dither_info.hpp
  src/core/dither_none.hpp
  src/core/dither_ordered.hpp
  src/core/frame_state.cpp
  src/core/frame_state.hpp
  src/core/grain.cpp
  src/core/grain.hpp
  src/core/kernel.hpp
  src/core/kernel_scalar.cpp
  src/core/kernel_scalar-inl.hpp
  src/core/kernel_scalar_dispatch.hpp
  src/core/kernel_scalar_inst.cpp
  src/core/kernel_hwy.cpp
  src/core/kernel_hwy-inl.hpp
  src/core/kernel_hwy_block-inl.hpp
  src/core/kernel_hwy_dither-inl.hpp
  src/core/kernel_hwy_modes-inl.hpp
  src/core/kernel_hwy_reference-inl.hpp
  src/core/kernel_hwy_store-inl.hpp
  src/core/math.hpp
  src/core/output_16bit.hpp
  src/core/parameters.cpp
  src/core/parameters.hpp
  src/core/pixel_proc.hpp
  src/core/pixel_proc_common.hpp
  src/core/plane.cpp
  src/core/plane.hpp
  src/core/processor.cpp
  src/core/processor.hpp
  src/core/random.cpp
  src/core/random.hpp
  src/f3kdb.h
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
  foreach(backend IN ITEMS purec highway)
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
