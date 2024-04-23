#! /usr/bin/bash

mkdir -p build
cd build
cmake -DCMAKE_C_FLAGS=-pg -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg ..
cmake --build .

