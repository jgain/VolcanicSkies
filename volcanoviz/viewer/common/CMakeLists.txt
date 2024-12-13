set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(COMMON_SOURCES
    initialize.cpp initialize.h
    mathutils.cpp mathutils.h
    progress.cpp progress.h
    region.cpp region.h
    stats.cpp stats.h
    timer.cpp timer.h
    cnpy.cpp cnpy.h
    )

if (BUILD_SOURCE2CPP)
    set(KERNELS
        ../terviewer/shaders/basic.vert
        ../terviewer/shaders/basic.frag
        ../terviewer/shaders/genNormal.vert
        ../terviewer/shaders/genNormal.frag
        ../terviewer/shaders/simple.vert
        ../terviewer/shaders/simple.frag
        ../terviewer/shaders/phong.vert
        ../terviewer/shaders/phong.frag
        ../terviewer/shaders/phongRS.vert
        ../terviewer/shaders/phongRS.frag
        ../terviewer/shaders/phongRSmanip.vert
        ../terviewer/shaders/phongRSmanip.frag
        ../terviewer/shaders/rad_scaling_pass1.vert
        ../terviewer/shaders/rad_scaling_pass1.frag
        ../terviewer/shaders/rad_scaling_pass2.vert
        ../terviewer/shaders/rad_scaling_pass2.frag
        ../terviewer/shaders/sun.vert
        ../terviewer/shaders/sun.frag
        ../terviewer/shaders/canopy.vert
        ../terviewer/shaders/canopy.frag)
    add_custom_command(
        OUTPUT source2cpp.cpp
        COMMAND ${PYTHON_EXECUTABLE} source2cpp.py ${KERNELS} ${CMAKE_CURRENT_BINARY_DIR}/source2cpp.cpp
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        MAIN_DEPENDENCY source2cpp.py
        DEPENDS ${KERNELS}
        COMMENT "Creating C++ source for shaders and kernels" VERBATIM)
    list(APPEND COMMON_SOURCES source2cpp.cpp)
endif()

add_library(common ${COMMON_SOURCES})
target_link_libraries(common ${Boost_SERIALIZATION_LIBRARY})
