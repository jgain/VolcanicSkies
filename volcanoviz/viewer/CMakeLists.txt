# This is needed for automoc
project(TerViewer)
cmake_minimum_required(VERSION 3.23)
include(CheckLibraryExists)
#set(CMAKE_PREFIX_PATH $ENV{HOME}/Qt/5.15.2/gcc_64)
#set(CMAKE_PREFIX_PATH $(brew --prefix qt5)/lib/cmake/Qt5)
#set(CMAKE_PREFIX_PATH $(brew --prefix qt5)/lib/cmake/Qt5Widgets)
#set(CMAKE_PREFIX_PATH $(brew --prefix qt5)/lib/cmake)
set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt@5/lib/cmake/Qt5")
#set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt@5/lib/cmake/Qt5")
set(Qt5Widgets_DIR $/opt/homebrew/opt/qt@5/lib/cmake/Qt5Widgets)

option(SMOKE_TEST "run a smoke test on every build" 1)
option(COVERAGE "compile for use with gcov" 0)
option(ASAN "compile with the address sanitiser" 0)
option(TSAN "compile with the thread sanitiser" 0)
option(SYNTHESIS_STATS "collect extra statistics about synthesis" 0)
enable_testing()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

### Check for dependencies
set(Boost_USE_MULTITHREADED ON)
add_definitions(-DBOOST_FILESYSTEM_VERSION=3)
find_package(Boost REQUIRED COMPONENTS filesystem program_options serialization system)
find_package(PkgConfig REQUIRED)
find_package(PythonInterp)
#find_package(Qt5Widgets 5 REQUIRED)
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
find_package(Qt5OpenGL 5 REQUIRED)
find_package(OpenGL)
find_package(ZLIB REQUIRED)
#find_package(GLUT)
pkg_check_modules(GLEW glew)

### Set include and link directories and definitions
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/khronos_headers
    ${Boost_INCLUDE_DIRS}
)
link_directories(
    ${PROJECT_SOURCE_DIR}/common
    ${PROJECT_SOURCE_DIR}/clh
    ${PROJECT_SOURCE_DIR}/terviewer
    ${Boost_LIB_DIRS}
)

if (Qt5Widgets_FOUND)
    include_directories(${Qt5Widgets_INCLUDE_DIRS})
    add_definitions(${Qt5Widgets_DEFINITIONS})
endif()
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_directories(${GLEW_LIBRARY_DIRS})
endif()
if (Qt5OpenGL_FOUND)
    include_directories(${Qt5OpenGL_INCLUDE_DIRS})
    add_definitions(${Qt5OpenGL_DEFINITIONS})
endif()
if (OPENGL_FOUND)
    include_directories(${OPENGL_INCLUDE_DIR})
endif()
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARY_DIRS})
endif()

### Check which components to build
if (PYTHONINTERP_FOUND AND OPENGL_FOUND)
    set(BUILD_SOURCE2CPP TRUE)
endif()
if (OPENCL_FOUND AND PYTHONINTERP_FOUND)
    set(BUILD_SYNTHESIZER TRUE)
    set(BUILD_CLH TRUE)
endif()
if (Qt5Widgets_FOUND
        AND Qt5OpenGL_FOUND
        AND OPENGL_FOUND
        AND GLEW_FOUND
        AND BUILD_SYNTHESIZER)
    set(BUILD_GUI TRUE)
endif()

### Set flags specific to the compiler
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
    # Eigen gives lots of warnings for this
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-enum-compare")
    # On Apple the framework headers are found ahead of the khronos_headers
    # directory, and they give deprecation warnings about OpenCL 1.1 functions.
    if (APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    endif()
    #if (NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_CXX_LINKER_FLAGS_RELEASE} -O4")
    #endif()
    # Link time optimization gives much better performance
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -flto")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -D_GLIBCXX_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} -flto")
endif()
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wno-unused-local-typedefs")
    if (${COVERAGE})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    endif()
    if (${ASAN})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libasan")
    endif()
    if (${TSAN})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -fpie -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pie -static-libtsan")
    endif()
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Eigen generates lots of warnings due to the register keyword being deprecated
    add_definitions(-Dregister=)
endif()
if (SYNTHESIS_STATS)
    add_definitions(-DSYNTHESIS_STATS)
endif()
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DUTS_DEBUG_CONTAINERS")


add_subdirectory(common)
add_subdirectory(clh)
add_subdirectory(terviewer)
