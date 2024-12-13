cmake_minimum_required(VERSION 3.5)

project(VolcanoWeather LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS true)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native -fopenmp -fPIC")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -mcpu=apple-m1 -Xpreprocessor -fopenmp -fPIC -I /usr/local/include/eigen3/ -I /opt/homebrew/opt/libomp/include")
#-I /opt/homebrew/Cellar/openvdb/10.0.1_1/include")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/OpenVDB")
#list(APPEND CMAKE_MODULE_PATH "/opt/homebrew/Cellar/openvdb/10.0.1_1/lib/cmake/OpenVDB")

set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY omp)

#set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lomp")
set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lOpenVDB")
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenVDB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Qt5 COMPONENTS Core Widgets 3DCore 3DRender 3DExtras REQUIRED)

include_directories(lib/)
include_directories(src/Engines)
include_directories(src/main)
include_directories(src/Models)
include_directories(src/Plume)
include_directories(src/PWMDataStructure)
include_directories(src/Utils)

link_directories("/opt/homebrew/opt/libomp/lib")
#link_directories("/opt/homebrew/Cellar/libomp/16.0.5/lib")

FILE(GLOB_RECURSE HEADER_FILES lib/*.h lib/*.hpp src/*.h src/*.hpp)

# QtCreator supports the following variables for Android, which are identical to qmake Android variables.
# Check http://doc.qt.io/qt-5/deployment-android.html for more information.
# They need to be set before the find_package(Qt5 ...) call.

#if(ANDROID)
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
#    if (ANDROID_ABI STREQUAL "armeabi-v7a")
#        set(ANDROID_EXTRA_LIBS
#            ${CMAKE_CURRENT_SOURCE_DIR}/path/to/libcrypto.so
#            ${CMAKE_CURRENT_SOURCE_DIR}/path/to/libssl.so)
#    endif()
#endif()

if(ANDROID)
  add_library(VolcanoWeather SHARED
      lib/json.hpp
      lib/OBJ_Loader.h
      src/Engines/abstractEngine.h
      src/Engines/advectionEngine.h
      src/Engines/atmosphericHeatingEngine.h
      src/Engines/heatDissipationEngine.h
      src/Engines/moistureUptakeEngine.h
      src/Engines/obstaclesEngine.h
      src/Engines/phaseTransitionEngine.h
      src/Engines/precipitationEngine.h
      src/Engines/pressureEngine.h
      src/Engines/sunEngine.h
      src/Engines/terrainHeatingEngine.h
      src/Engines/verticalCouplingEngine.h
      src/Engines/worldEngine.h

      src/main/main.cpp
      src/main/mainwindow.h
      src/main/mainwindow.cpp
      src/main/mainwindow.ui
  #    src/main/volcanoWeather.h
  #    src/main/volcanoWeather.cpp
      src/main/displayWindow.h
      src/main/displayWindow.cpp

      src/Models/airLayer.h
      src/Models/convectionLayer.h
      src/Models/planet.h
      src/Models/sun.h
      src/Models/terrain.h
      src/Models/world.h

      src/Plume/smoke.hpp
      src/Plume/skirt.hpp
      src/Plume/timer.hpp
      src/Plume/smokeLayer.hpp
      src/Plume/terrain_structure.hpp

      src/PWMDataStructure/AbstractPWMDataStructure.h
      src/PWMDataStructure/Coordinate.h
      src/PWMDataStructure/position.h
      src/PWMDataStructure/square2DArray.h

      src/Utils/atmofuncs.h
      src/Utils/cloudUtils.h
      src/Utils/mathUtils.h
      src/Utils/objUtils.h

      src/Utils/settings.h
      src/Utils/vdbExporter.h
      src/Utils/vizUtils.h

      src/Engines/abstractEngine.cpp
      src/PWMDataStructure/Coordinate.cpp
      src/Models/planet.cpp
      src/Utils/polynomials.cpp
      src/Utils/settings.cpp
      src/Plume/smoke.cpp
      src/Plume/skirt.cpp
      src/Plume/timer.cpp
      src/Plume/terrain_structure.cpp
  )
else()
    add_executable(VolcanoWeather
    lib/json.hpp
    lib/OBJ_Loader.h
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/atmosphericHeatingEngine.h
    src/Engines/heatDissipationEngine.h
    src/Engines/moistureUptakeEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/phaseTransitionEngine.h
    src/Engines/precipitationEngine.h
    src/Engines/pressureEngine.h
    src/Engines/sunEngine.h
    src/Engines/terrainHeatingEngine.h
    src/Engines/verticalCouplingEngine.h
    src/Engines/worldEngine.h

    src/main/main.cpp
    src/main/mainwindow.h
    src/main/mainwindow.cpp
    src/main/mainwindow.ui
#    src/main/volcanoWeather.h
#    src/main/volcanoWeather.cpp
    src/main/displayWindow.h
    src/main/displayWindow.cpp

    src/Models/airLayer.h
    src/Models/convectionLayer.h
    src/Models/planet.h
    src/Models/sun.h
    src/Models/terrain.h
    src/Models/world.h

    src/Plume/smoke.hpp
    src/Plume/skirt.hpp
    src/Plume/timer.hpp
    src/Plume/smokeLayer.hpp
    src/Plume/terrain_structure.hpp

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/position.h
    src/PWMDataStructure/square2DArray.h

    src/Utils/atmofuncs.h
    src/Utils/cloudUtils.h
    src/Utils/mathUtils.h
    src/Utils/objUtils.h

    src/Utils/settings.h
    src/Utils/vizUtils.h

    src/Engines/abstractEngine.cpp
    src/PWMDataStructure/Coordinate.cpp
    src/Models/planet.cpp
    src/Utils/settings.cpp
    src/Plume/smoke.cpp
    src/Plume/skirt.cpp
    src/Plume/timer.cpp
    src/Plume/terrain_structure.cpp
  )
endif()


#add_executable(CoordinateTest ${HEADER_FILES}
#    src/Tests/CoordinateTest.cpp
#    src/PWMDataStructure/Coordinate.cpp
#)

#add_executable(PositionTest ${HEADER_FILES}
#    src/Tests/PositionTest.cpp
#    src/PWMDataStructure/position.cpp
#)

#add_executable(Square2DArrayTest ${HEADER_FILES}
#    src/Tests/Square2DArrayTest.cpp
#    src/PWMDataStructure/Coordinate.cpp
#)

#add_executable(AirLayerTest ${HEADER_FILES}
#    src/Tests/AirLayerTest.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Models/planet.cpp
#)

#add_executable(AdvectionEngineTest
#    src/Engines/abstractEngine.h
#    src/Engines/advectionEngine.h
#    src/Engines/pressureEngine.h

#    src/Models/airLayer.h
#    src/Models/planet.h

#    src/PWMDataStructure/AbstractPWMDataStructure.h
#    src/PWMDataStructure/Coordinate.h
#    src/PWMDataStructure/square2DArray.h
#    src/PWMDataStructure/flatStaggeredGrid.h

#    src/Utils/atmofuncs.h
#    src/Utils/mathUtils.h

#    src/Utils/vdbExporter.h
#    src/Utils/vizUtils.h

#    src/Tests/AdvectionEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#    src/Models/planet.cpp
#    src/PWMDataStructure/Coordinate.cpp
#)
#target_link_libraries(AdvectionEngineTest "/opt/homebrew/opt/libomp/lib/libomp.a")
#target_link_libraries(AdvectionEngineTest ${OpenCV_LIBS})

add_executable(FlatMACTest
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/pressureEngine.h

    src/Models/airLayer.h
    src/Models/planet.h

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/square2DArray.h
    src/PWMDataStructure/flatStaggeredGrid.h

    src/Utils/atmofuncs.h
    src/Utils/mathUtils.h

    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Tests/FlatMACTest.cpp
    src/Engines/abstractEngine.cpp
    src/Models/planet.cpp
    src/PWMDataStructure/Coordinate.cpp
)
#target_link_libraries(FlatMACTest "/opt/homebrew/opt/libomp/lib/libomp.a")
target_link_libraries(FlatMACTest ${OpenCV_LIBS})

add_executable(FlatMACUlyssesComparison
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/pressureEngine.h

    src/Models/airLayer.h
    src/Models/planet.h

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/square2DArray.h
    src/PWMDataStructure/flatStaggeredGrid.h

    src/Utils/atmofuncs.h
    src/Utils/mathUtils.h

    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Tests/FlatMACUlyssesComparison.cpp
    src/Engines/abstractEngine.cpp
    src/Models/planet.cpp
    src/PWMDataStructure/Coordinate.cpp
)
#target_link_libraries(FlatMACUlyssesComparison "/opt/homebrew/opt/libomp/lib/libomp.a")
target_link_libraries(FlatMACUlyssesComparison ${OpenCV_LIBS})

#add_executable(PlanetTest ${HEADER_FILES}
#    src/Tests/PlanetTest.cpp
#    src/Models/planet.cpp
#)

#add_executable(TerrainTest ${HEADER_FILES}
#    src/Tests/TerrainTest.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Models/planet.cpp
#)

#add_executable(ConvectionLayerTest ${HEADER_FILES}
#    src/Engines/abstractEngine.h
#    src/Engines/advectionEngine.h
#    src/Engines/obstaclesEngine.h
#    src/Engines/verticalCouplingEngine.h

#    src/Models/airLayer.h
#    src/Models/planet.h

#    src/PWMDataStructure/AbstractPWMDataStructure.h
#    src/PWMDataStructure/Coordinate.h
#    src/PWMDataStructure/square2DArray.h
#    src/PWMDataStructure/flatStaggeredGrid.h

#    src/Utils/atmofuncs.h
#    src/Utils/mathUtils.h

#    src/Utils/vdbExporter.h
#    src/Utils/vizUtils.h

#    src/Engines/abstractEngine.cpp
#    src/Tests/ConvectionLayerTest.cpp
#    src/Models/planet.cpp
#    src/PWMDataStructure/Coordinate.cpp
#)

add_executable(VerticalCouplingEngineTest
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/pressureEngine.h
    src/Engines/verticalCouplingEngine.h

    src/Models/airLayer.h
    src/Models/planet.h

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/square2DArray.h
    src/PWMDataStructure/flatStaggeredGrid.h

    src/Utils/atmofuncs.h
    src/Utils/mathUtils.h

    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Tests/VerticalCouplingEngineTest.cpp
    src/Engines/abstractEngine.cpp
    src/Models/planet.cpp
    src/PWMDataStructure/Coordinate.cpp
)
#target_link_libraries(VerticalCouplingEngineTest "/opt/homebrew/opt/libomp/lib/libomp.a")
target_link_libraries(VerticalCouplingEngineTest ${OpenCV_LIBS})

#add_executable(WorldTest ${HEADER_FILES}
#    src/Tests/WorldTest.cpp
#    src/Models/planet.cpp
#    src/Plume/terrain_structure.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Utils/settings.cpp
#)

#add_executable(HeatDissipationEngineTest ${HEADER_FILES}
#    src/Tests/HeatDissipationEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#    src/Models/planet.cpp
#    src/PWMDataStructure/Coordinate.cpp
#)

#add_executable(PhaseTransitionEngineTest ${HEADER_FILES}
#    src/Tests/PhaseTransitionEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#    src/Models/planet.cpp
#    src/PWMDataStructure/Coordinate.cpp
#)

#add_executable(SunTest ${HEADER_FILES}
#    src/Tests/SunTest.cpp
#)

#add_executable(SunEngineTest ${HEADER_FILES}
#    src/Tests/SunEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#)

#add_executable(TerrainHeatingEngineTest ${HEADER_FILES}
#    src/Tests/TerrainHeatingEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Models/planet.cpp
#)

#add_executable(MoistureUptakeEngineTest ${HEADER_FILES}
#    src/Tests/MoistureUptakeEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Models/planet.cpp
#)

#add_executable(AtmosphericHeatingEngineTest ${HEADER_FILES}
#    src/Tests/AtmosphericHeatingEngineTest.cpp
#    src/Engines/abstractEngine.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Models/planet.cpp
#)

add_executable(WorldEngineTest
    lib/json.hpp
    lib/OBJ_Loader.h
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/atmosphericHeatingEngine.h
    src/Engines/heatDissipationEngine.h
    src/Engines/moistureUptakeEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/phaseTransitionEngine.h
    src/Engines/precipitationEngine.h
    src/Engines/pressureEngine.h
    src/Engines/sunEngine.h
    src/Engines/terrainHeatingEngine.h
    src/Engines/verticalCouplingEngine.h
    src/Engines/worldEngine.h

    src/Models/airLayer.h
    src/Models/convectionLayer.h
    src/Models/planet.h
    src/Models/sun.h
    src/Models/terrain.h
    src/Models/world.h

    src/Plume/terrain_structure.hpp

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/position.h
    src/PWMDataStructure/square2DArray.h
    src/PWMDataStructure/flatStaggeredGrid.h

    src/Utils/atmofuncs.h
    src/Utils/cloudType.h
    src/Utils/cloudUtils.h
    src/Utils/mathUtils.h
    src/Utils/objUtils.h

    src/Utils/settings.h
    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Tests/WorldEngineTest.cpp
    src/Engines/abstractEngine.cpp
    src/PWMDataStructure/Coordinate.cpp
    src/Models/planet.cpp
    src/Plume/terrain_structure.cpp
)
#target_link_libraries(WorldEngineTest "/opt/homebrew/opt/libomp/lib/libomp.a")
target_link_libraries(WorldEngineTest ${OpenCV_LIBS})

add_executable(MergeTestbed
    lib/json.hpp
    lib/OBJ_Loader.h
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/atmosphericHeatingEngine.h
    src/Engines/heatDissipationEngine.h
    src/Engines/moistureUptakeEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/phaseTransitionEngine.h
    src/Engines/precipitationEngine.h
    src/Engines/pressureEngine.h
    src/Engines/sunEngine.h
    src/Engines/terrainHeatingEngine.h
    src/Engines/verticalCouplingEngine.h
    src/Engines/worldEngine.h

    src/Models/airLayer.h
    src/Models/convectionLayer.h
    src/Models/planet.h
    src/Models/sun.h
    src/Models/terrain.h
    src/Models/world.h

    src/Plume/smoke.hpp
    src/Plume/skirt.hpp
    src/Plume/timer.hpp
    src/Plume/smokeLayer.hpp
    src/Plume/terrain_structure.hpp

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/position.h
    src/PWMDataStructure/square2DArray.h
    src/PWMDataStructure/flatStaggeredGrid.h

    src/Utils/atmofuncs.h
    src/Utils/cloudType.h
    src/Utils/cloudUtils.h
    src/Utils/mathUtils.h
    src/Utils/objUtils.h

    src/Utils/settings.h
    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Engines/abstractEngine.cpp
    src/PWMDataStructure/Coordinate.cpp
    src/Models/planet.cpp
    src/Utils/settings.cpp
    src/Plume/smoke.cpp
    src/Plume/skirt.cpp
    src/Plume/timer.cpp
    src/Plume/terrain_structure.cpp
    src/Tests/MergeTest.cpp
)
#target_link_libraries(MergeTestbed "/opt/homebrew/opt/libomp/lib/libomp.a")

#add_executable(CloudClassificationTest ${HEADER_FILES}
#    src/Tests/CloudClassificationTest.cpp
#    src/PWMDataStructure/Coordinate.cpp
#    src/Models/planet.cpp
#    src/Plume/terrain_structure.cpp
#)

#add_executable(SettingsTest ${HEADER_FILES}
#    src/Tests/SettingsTest.cpp
#    src/Models/planet.cpp
#    src/Utils/settings.cpp
#)

add_executable(ObstaclesEngineTest
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/pressureEngine.h
    src/Engines/verticalCouplingEngine.h
    src/Models/airLayer.h
    src/Models/convectionLayer.h
    src/Models/planet.h
    src/Models/terrain.h

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/square2DArray.h
    src/PWMDataStructure/flatStaggeredGrid.h

    src/Utils/atmofuncs.h
    src/Utils/cloudType.h
    src/Utils/cloudUtils.h
    src/Utils/mathUtils.h
    src/Utils/objUtils.h

    src/Utils/settings.h
    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Tests/ObstaclesEngineTest.cpp
    src/Engines/abstractEngine.cpp
    src/PWMDataStructure/Coordinate.cpp
    src/Models/planet.cpp
    src/Plume/terrain_structure.cpp
)
#target_link_libraries(ObstaclesEngineTest "/opt/homebrew/opt/libomp/lib/libomp.a")
target_link_libraries(ObstaclesEngineTest ${OpenCV_LIBS})

add_executable(PrecipitationEngineTest
    src/Engines/abstractEngine.h
    src/Engines/advectionEngine.h
    src/Engines/obstaclesEngine.h
    src/Engines/precipitationEngine.h
    src/Engines/pressureEngine.h
    src/Engines/verticalCouplingEngine.h
    src/Models/airLayer.h
    src/Models/convectionLayer.h
    src/Models/planet.h
    src/Models/terrain.h

    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/square2DArray.h

    src/Utils/atmofuncs.h
    src/Utils/cloudType.h
    src/Utils/cloudUtils.h
    src/Utils/mathUtils.h
    src/Utils/objUtils.h

    src/Utils/settings.h
    src/Utils/vdbExporter.h
    src/Utils/vizUtils.h

    src/Tests/PrecipitationEngineTest.cpp
    src/Engines/abstractEngine.cpp
    src/PWMDataStructure/Coordinate.cpp
    src/Models/planet.cpp
    src/Plume/terrain_structure.cpp
)
#target_link_libraries(PrecipitationEngineTest "/opt/homebrew/opt/libomp/lib/libomp.a")

add_executable(SmoothingTest
    src/PWMDataStructure/AbstractPWMDataStructure.h
    src/PWMDataStructure/Coordinate.h
    src/PWMDataStructure/Coordinate.cpp
    src/PWMDataStructure/flatStaggeredGrid.h
    src/Tests/SmoothingTest.cpp
)

target_link_libraries(MergeTestbed OpenVDB::openvdb Eigen3::Eigen ${OpenCV_LIBS})
target_link_libraries(PrecipitationEngineTest OpenVDB::openvdb Eigen3::Eigen ${OpenCV_LIBS})
#target_link_libraries(CloudClassificationTest OpenVDB::openvdb Eigen3::Eigen)
target_link_libraries(VolcanoWeather PRIVATE Qt5::Widgets Qt5::3DCore Qt5::3DRender Qt5::3DExtras Eigen3::Eigen)
