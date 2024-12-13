#!/bin/sh

rm -r ./xcode
mkdir xcode
rm CMakeCache.txt
cd xcode
#cmake -DQt5_DIR=$(brew --prefix qt5)/lib/cmake/Qt5 ..
#cmake -DQt5Widgets_DIR=$/opt/homebrew/opt/qt@5/lib/cmake/Qt5Widgets ..
#cmake -DGLUT_ROOT=$/opt/homebrew/opt/freeglut ..
#cmake -DCMAKE_CXX_COMPILER=g++ ..
#cmake -DCMAKE_BUILD_TYPE=Release ..
cmake -G Xcode ..
#make
