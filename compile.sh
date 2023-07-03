#!/bin/bash

rm -rf build
rm -rf lib
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install -j 8
cd ..
