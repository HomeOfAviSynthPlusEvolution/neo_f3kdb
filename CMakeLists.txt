cmake_minimum_required (VERSION 3.1)
project (neo_f3kdb)
file(GLOB AVS_FILES ${PROJECT_SOURCE_DIR}/src/avisynth/*.cpp)
file(GLOB VPY_FILES ${PROJECT_SOURCE_DIR}/src/vapoursynth/*.cpp)
file(GLOB CORE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
find_package(Git)
add_library(neo-f3kdb SHARED ${VPY_FILES} ${CORE_FILES} version.rc)
if(MSVC)
  target_sources(neo-f3kdb PRIVATE ${AVS_FILES})
endif(MSVC)
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
add_compile_definitions(FLASH3KYUU_DEBAND_EXPORTS)
execute_process(COMMAND ${GIT_EXECUTABLE} describe --first-parent --tags --always OUTPUT_VARIABLE GIT_REPO_VERSION)
string(REGEX REPLACE "(r[0-9]+).*" "\\1" VERSION ${GIT_REPO_VERSION})
configure_file (
  "${PROJECT_SOURCE_DIR}/version.rc.in"
  "${PROJECT_SOURCE_DIR}/version.rc"
)
include_directories(include)
include_directories(${PROJECT_SOURCE_DIR}/src/avisynth)
include_directories(${PROJECT_SOURCE_DIR}/src/vapoursynth)
add_custom_command(
  TARGET neo-f3kdb POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:neo-f3kdb> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:neo-f3kdb>"
)
