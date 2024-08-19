cmake_minimum_required(VERSION 3.15)
project(neo_f3kdb CXX)

file(GLOB CODE "src/*.cpp")
file(GLOB SSE41_CODE_IMPL "src/*sse4.cpp")
add_library(neo-f3kdb SHARED main.cpp src/version.rc ${CODE} ${SSE41_CODE_IMPL})
set_property(TARGET neo-f3kdb PROPERTY CXX_STANDARD 17)
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
    include_directories(${AVISYNTH_INCLUDE_DIRS})
  else()
    include_directories(include/avisynth)
  endif()

  pkg_check_modules(VAPOURSYNTH vapoursynth)
  if(VAPOURSYNTH_FOUND)
    include_directories(${VAPOURSYNTH_INCLUDE_DIRS})
  else()
    include_directories(include/vapoursynth)
  endif()
else()
  include_directories(include/avisynth)
  include_directories(include/vapoursynth)
endif()

include_directories(.)
include_directories(include/dualsynth)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
    set_source_files_properties(${SSE41_CODE_IMPL} PROPERTIES COMPILE_FLAGS "-march=penryn")
  else()
    set_source_files_properties(${SSE41_CODE_IMPL} PROPERTIES COMPILE_FLAGS "/arch:SSE2")
  endif()

  if(CMAKE_GENERATOR_TOOLSET MATCHES "v[0-9]*_xp")
    target_compile_definitions(neo-f3kdb PRIVATE WINVER=0x502 _WIN32_WINNT=0x502)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:threadSafeInit-")
  else()
    target_compile_options(neo-f3kdb PRIVATE /O2)
  endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(${SSE41_CODE_IMPL} PROPERTIES COMPILE_FLAGS "/arch:SSE2")
  target_link_libraries(neo-f3kdb libmmds)
  target_compile_options(neo-f3kdb PRIVATE /O3)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  set_source_files_properties(${SSE41_CODE_IMPL} PROPERTIES COMPILE_FLAGS "-msse4.1")
  target_link_libraries(neo-f3kdb libmmds)
  target_compile_options(neo-f3kdb PRIVATE -O3)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    target_compile_options(neo-f3kdb PRIVATE -O3 -funroll-loops)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
    target_compile_options(neo-f3kdb PRIVATE -O3 -funroll-loops)
    set_source_files_properties(${SSE41_CODE_IMPL} PROPERTIES COMPILE_FLAGS "-msse4.1")
  else()
    target_compile_options(neo-f3kdb PRIVATE -O3 -funroll-loops)
  endif()

else()
  target_compile_options(neo-f3kdb PRIVATE -O3)
endif()

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(execution HAS_EXECUTION)
if(HAS_EXECUTION)
  add_definitions(-DHAS_EXECUTION)
endif()
if(ENABLE_PAR AND HAS_EXECUTION)
  add_definitions(-DENABLE_PAR)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(neo-f3kdb tbb)
  endif()
endif()

add_custom_command(
  TARGET neo-f3kdb POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:neo-f3kdb> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:neo-f3kdb>"
)
