mkdir -p build/gcc
pushd build/gcc
cmake -G "MSYS Makefiles" -DCMAKE_CXX_FLAGS=-msse4.1 -D_DIR=gcc ../../
popd
cmake --build build/gcc
