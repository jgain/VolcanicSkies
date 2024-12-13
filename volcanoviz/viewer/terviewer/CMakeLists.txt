set (BUILD_GUI 1)
if (BUILD_GUI)
    
   set(CMAKE_AUTOMOC TRUE)
   set(CMAKE_INCLUDE_CURRENT_DIR TRUE)  # MOC files generated in the binary dir
   set(DATA_IMPORT_DIR ${PROJECT_SOURCE_DIR}/../data_importer)
   set(BASE_ALL_DIR ${PROJECT_SOURCE_DIR}/..)

   set(GUI_SOURCES
       glwidget.cpp glwidget.h
       timer.cpp timer.h
       descriptor.cpp descriptor.h
       dice_roller.cpp dice_roller.h
       mesh.cpp mesh.h
       noisefield.cpp noisefield.h
       palette.h palette.cpp
       stroke.cpp stroke.h
       terrain.cpp terrain.h
       volcano.cpp volcano.h
       shape.cpp shape.h
       typemap.cpp typemap.h
       vecpnt.cpp vecpnt.h
       view.cpp view.h
       waterfill.cpp waterfill.h
       window.cpp window.h
       shaderProgram.cpp shaderProgram.h
       trenderer.cpp trenderer.h
       ${DATA_IMPORT_DIR}/data_importer.cpp ${DATA_IMPORT_DIR}/data_importer.h
       ${BASE_ALL_DIR}/common/basic_types.h
)

   #set(CMAKE_BUILD_TYPE Debug)
   #set(CMAKE_BUILD_TYPE Release)
    add_library(vgui ${GUI_SOURCES})
    target_include_directories(vgui PUBLIC ${PROJECT_SOURCE_DIR} ${BASE_ALL_DIR})
    target_link_libraries(vgui common
        ${GLEW_LIBRARIES} ${GLUT_LIBRARIES} ${Qt5Widgets_LIBRARIES} ${Qt5OpenGL_LIBRARIES}
        ${OPENGL_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${ZLIB_LIBRARIES})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    add_executable(terviewer main.cpp)
    target_include_directories(terviewer PRIVATE ${PROJECT_SOURCE_DIR} ${BASE_ALL_DIR})
    target_link_libraries(terviewer vgui sqlite3) # stdc++fs)
endif()
