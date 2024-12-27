MAC INSTALLATION:
-----------------

Install HomeBrew if you haven't already (https://brew.sh).

Intel processors:

with M1 Processors:
brew install gcc 

Create a symbolic link from gcc-11 and g++-11, or whichever version you have just installed, to gcc and g++ to override the default use of clang

Something like:
cd /opt/homebrew/bin/
ln -s gcc-11 gcc
ln -s g++-11 g++

brew install qt5 mesa-glu freeglut glew zlib glm

It is fine to use gcc-11 

Uncomment the qt5 path setting in herdvid/buildherdapple.sh

UBUNTU INSTALLATION:
--------------------
These are requirements for Ubuntu. In each case, I've listed where to get the software 
for Ubuntu 18.10. For newer versions of Ubuntu you might not need all the PPAs.

CMake 2.8.7+, make, automake 1.9, pkg-config, GLUT, QT5, Boost, at least 1.49 (earlier versions don't play nice with C++11)::

sudo apt-get install cmake make automake1.9 pkg-config freeglut3-dev libglm-dev qtbase5-dev libboost-all-dev libglew-dev libcppunit-dev libomp-dev qtcreator sqlite3 sqlitebrowser libzmq3-dev

To install a particular version of g++ (in this case v7, note that version 9 does not work):

sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt install gcc-8 g++-8
then make sure to set the correct version in the build script (buildsim.sh)

ADDED PACKAGES:

sudo apt-get install libgl1-mesa-dev libqtcore4 libqtgui4 libqt4-opengl libqt4-opengl-dev //line used in 3rd year to fix opengl stuff 

sudo apt-get install libxmu-dev libxi-dev


ADDITIONAL REQUIREMENTS:

Modify the `basedir` variable in viewer/CMakeLists.txt to contain the location of this directory on your machine.

packages/libraries:
+ SDL2		(install: sudo apt install libsdl2-dev)
+ Assimp	(install: sudo apt install libassimp-dev)
+ pthread, X11, libpng, libz (although I think the first two are usually installed by default)

Compiling and Executing
-----------------------

#### Ubuntu & Apple

There is a build script that you can run from ./herdvid: buildherd.sh

Alternatively, once all the requirements are running, create a subdirectory called build (or
anything starting with build - you can have multiple build directories), switch
into it, and run

cmake <options> ..

Some useful options
-DCMAKE_CXX_COMPILER=g++-4.8          (to force a specific compiler)
-DCMAKE_BUILD_TYPE=Debug|Release|RelWithDebInfo  (to select build type, only choose one option here)
-DCMAKE_CXX_FLAGS_DEBUG="-gdwarf-3"   (helps older versions of GDB interpret debug symbols)
-DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O3 -DNDEBUG -g -gdwarf-3" (similar)
-DCMAKE_CXX_FLAGS="-fPIC" (if the compilation complains about position independent code)

Then run make to build. cmake options are sticky.  System must be run from the build directory because there are some relative paths.

Command line execution is: ./terviewer/terviewer <data directory> (where <data directory> should contain the terrain elevation file). Data samples are available in the directory /data

A good test execution is: ./terviewer/terviewer ../../data/S4500-4500-1024/

GUI
---
Once a model has been loaded, it can be rotated by holding down the right mouse button and moving the cursor. Zooming is with the mouse wheel. Double click with the right mouse button to change the focal point on the terrain.
