cmake_minimum_required(VERSION 3.15)
project(neo_f3kdb CXX)
file(GLOB CODE "src/*.cpp")
add_library(neo-f3kdb SHARED main.cpp src/version.rc ${CODE} ${CODE_IMPL})
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

if(NOT MSVC)
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

if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  target_compile_options(neo-f3kdb PRIVATE -Ofast -ftree-vectorize -fopenmp)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
  target_compile_options(neo-f3kdb PRIVATE -Ofast -msse4.2 -mavx -funroll-loops -fopenmp)
else()
  target_compile_options(neo-f3kdb PRIVATE -Ofast -funroll-loops -fopenmp)
endif()

add_custom_command(
  TARGET neo-f3kdb POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:neo-f3kdb> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:neo-f3kdb>"
)