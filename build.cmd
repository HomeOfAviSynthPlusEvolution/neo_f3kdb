@echo off

mkdir build\msvc-x86
pushd build\msvc-x86
cmake -DCMAKE_GENERATOR_PLATFORM=Win32 -D_DIR=msvc-x86 ..\..\
popd
mkdir build\msvc-x64
pushd build\msvc-x64
cmake -DCMAKE_GENERATOR_PLATFORM=x64 -D_DIR=msvc-x64 ..\..\
popd
cmake --build build\msvc-x86 --config Release
cmake --build build\msvc-x64 --config Release

mkdir build\icc-x86
pushd build\icc-x86
cmake -T"Intel C++ Compiler 19.0" -DCMAKE_GENERATOR_PLATFORM=Win32 -D_DIR=icc-x86 ..\..\
popd
mkdir build\icc-x64
pushd build\icc-x64
cmake -T"Intel C++ Compiler 19.0" -DCMAKE_GENERATOR_PLATFORM=x64 -D_DIR=icc-x64 ..\..\
popd
cmake --build build\icc-x86 --config Release
cmake --build build\icc-x64 --config Release
