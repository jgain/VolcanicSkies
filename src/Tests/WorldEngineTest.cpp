/*
 * Volcanic Skies
 * Copyright (C) 2024 P. Cilliers Pretorius, University of Cape Town
 *
 * This file is part of the Volcanic Skies project.
 *
 * Volcanic Skies is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License (GPL) as published by the Free Software 
 * Foundation, either version 2 of the License, or (at your discretion) any later version.
 *
 * Volcanic Skies is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with 
 * Volcanic Skies in the LICENSE file. If not, see <https://www.gnu.org/licenses/>.
 *
 * Additional information and disclaimers regarding liability and third-party 
 * components can be found in the NOTICE file included with this project.
 *
 */
#include "cloudUtils.h"
#include <filesystem>
#include <iostream>
#include <limits>
#include "vizUtils.h"
#include "worldEngine.h"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;
using namespace PWM::Engine;

terrain_structure constructTestTerrStruct(std::string file, float cellWidth = 400.0f, float sceneWidth = 15000.0f){
    terrain_structure res(file, cellWidth, sceneWidth);
    return res;
}

void testConstructors(){
	auto x = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj");
    x->init("../resources/PlanetEarth.json", ter);
	
	auto a = worldEngine(x, 60.f, true);
	if (a.getWorldModel() == x)
		std::cout << "\033[1;32mConstruction works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Construction does not perform as expected!\033[0m" << std::endl;
}

void testGetters(){
	auto x = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj");
    x->init("../resources/PlanetEarth.json", ter);

	auto a = worldEngine(x, 60.f, true);
	if (a.getWorldModel() == x)
		std::cout << "\033[1;32mgetWorldModel() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getWorldModel() does not perform as expected!\033[0m" << std::endl;
}

void testMiscellany(int layers, int width, int height, int steps, valType xSize, valType ySize, valType topOfAtmo){
	/* auto x = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
	x->init("../resources/testResources/PlanetEarth.json", "../resources/testResources/test.pwea");
	
	auto a = worldEngine(x, 60.f, true);

	for (int i = 0; i < 1000; ++i)
		a.step();

	std::cout << "Time taken for 1000 steps of WorldEngine for a 6 x 6 grid: " << a.getRunTimePassed() << " seconds." << std::endl;
	*/
    auto currDir = std::filesystem::current_path();
	for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
		if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".png")
			std::filesystem::remove(dir_Entry.path());
    }

    currDir.remove_filename().remove_filename();
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir.concat("ImageOutput/WorldEngineTesting/ConstantWind");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir.remove_filename().concat("RandomWind");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

//	auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto z = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj", 400.f, xSize / 2);
//    auto ter = terrain_structure();
//    ter.testTerrainStructure(width);
//    y->maxLayers = layers;
//    y->atmoTop = topOfAtmo;
    z->maxLayers = layers;
    z->atmoTop = topOfAtmo;

//    y->init("../resources/PlanetEarth.json", ter, width, height, xSize, ySize);
//    y->addSun("../resources/SunSol.json");
    z->init("../resources/PlanetEarth.json", ter, width, height, xSize, ySize);
    z->addSun("../resources/SunSol.json");

    for (int i = 0; i < layers; ++i){
//        y->getAirLayer(i)->getVelocityTheta() = 8;
//        y->getAirLayer(i)->getVelocityPhi() = 8;
        z->getAirLayer(i)->randomInit();
//        y->getAirLayer(i)->getMoisture() = 0.05;
        z->getAirLayer(i)->getMoisture() = 0.01;
//        y->getAirLayer(i)->getCondensedWater() = 0;
        z->getAirLayer(i)->getCondensedWater() = 0;
        z->getAirLayer(i)->getTemperature() = 300 - 10 * i;
        z->getAirLayer(i)->getVelocityTheta() = 0;
        z->getAirLayer(i)->getVelocityPhi() = 0;
    }

//    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
//    int hotWidth = width / 10;
//    int startj = height / 2, startk = width / 2;
//    double highT = 340.0;
//    for (int j = -hotWidth; j < hotWidth; ++j){
//        for (int k = -hotWidth; k < hotWidth; ++k){
//            y->getTerrain()->setTemperature(ijToI(width, j + startj, k + startk), highT);
//            z->getTerrain()->setTemperature(ijToI(width, j + startj, k + startk), highT);
//        }
//    }

//    y->getAirLayer(3)->getCondensedWater() = 0.2f;
    z->getAirLayer(3)->getTemperature() = 245;
    z->getAirLayer(3)->getCondensedWater() = 0.12f;

//    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/WorldEngineTesting/PWMElevationMap" + std::to_string(y->getTerrain()->getElevation().getX()) + "x" + std::to_string(y->getTerrain()->getElevation().getY()) + ".ppm", y->getTerrain()->getElevation());
    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/WorldEngineTesting/PlumeElevationMap" + std::to_string(ter.field_size) + "x" + std::to_string(ter.field_size) + ".ppm", ter);
//    std::cout << "Lava lake is at " << y->getTerrain()->getElevation((height / 2) * width + width / 2) << " m above ASL." << std::endl;
    std::stringstream terTF;
    terTF << "../ImageOutput/WorldEngineTesting/Terrain_Init_Temperature.ppm";
//    PWM::Utils::writeTempImage(terTF.str(), y->getTerrain()->getTemperature());

//    auto b = worldEngine(y, 1.f, true);
    auto c = worldEngine(z, 1.f, true);

//	if (b.getWorldModel() == y)
//		std::cout << "\033[1;32mUsing a real world map works as expected.\033[0m" << std::endl;
//	else
//		std::cout << "\033[1;31mError! Using a real world map does not perform as expected!\033[0m" << std::endl;
	
    double previousStep = 0;

	auto manCloudGen = [](const std::shared_ptr<PWM::Model::airLayer<dsType, valType>>& layer, const size_t gap, const valType val){
		#pragma omp parallel for
			for (size_t i = 0; i < layer->getObstacles().size(); ++i){
				if ((i / gap) % 5 == 0)
					layer->setCondensedWater(i, val);
				else if ((i / gap) % 5 == 1)
					layer->setCondensedWater(i, val / 100);
			}
		#pragma omp barrier
	};
//	std::cerr << "Testing manually set clouds in cloud extraction" << std::endl;
//	manCloudGen(b.getWorldModel()->getAirLayer(0), 20, 4);
//	manCloudGen(b.getWorldModel()->getAirLayer(2), 40, 4);
//	manCloudGen(b.getWorldModel()->getAirLayer(4), 80, 4);
//	manCloudGen(b.getWorldModel()->getAirLayer(6), 300, 4);
//	manCloudGen(b.getWorldModel()->getAirLayer(b.getWorldModel()->getAirLayers().size() - 1), 125, 4);
//  PWM::Utils::extractClouds(b.getWorldModel(), "TestingOutput", 0.1);
//	PWM::Utils::basicCloudExtract(b.getWorldModel(), "TestingBasicOutput", 0.1);
//    PWM::Utils::writeAirImages("../ImageOutput/WorldEngineTesting", 0, b.getWorldModel());
//    b.printViz = false;
    c.printViz = false;
//    b.outputDirectory = "../ImageOutput/WorldEngineTesting/ConstantWind";
    c.outputDirectory = "../ImageOutput/WorldEngineTesting/RandomWind";
    for (int i = 0; i <= steps; ++i){
        if (i % 10 == 0){
//            PWM::Utils::writeAirImages("../ImageOutput/WorldEngineTesting/ConstantWind", i, b.getWorldModel());
//            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/WorldEngineTesting/ConstantWind", i, 0.1);
            PWM::Utils::writeAirImages("../ImageOutput/WorldEngineTesting/RandomWind", i, c.getWorldModel());
            PWM::Utils::basicCloudExtract(c.getWorldModel(), "../ImageOutput/WorldEngineTesting/RandomWind", i, 0.1);
        }
//        previousStep = b.getRunTimePassed();
//        b.step();
        c.step();
//        double currentStep = b.getRunTimePassed();
//        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;
        if (i == steps / 8)
            std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
        else if (i == steps / 4)
            std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 8))
            std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
        else if (i == steps / 2)
            std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
        else if (i == (5 * steps / 8))
            std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 4))
            std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
        else if (i == (7 * steps / 8))
            std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
	}
//    PWM::Utils::writeAirGifImages("ImageOutput/WorldEngineTesting", b.getWorldModel());
//    b.printTiming(true);
    c.printTiming(true);
}

int main(int argc, char** argv){
//    omp_set_num_threads(1);
    int layers = 8;
    int reps = 10;
    int width = 512;
    int height = 512;
    int steps = 600;
    switch (argc){
            case 1 :
                std::cout << "Configuration is available through program arguments in the form \"./WorldEngineTest [#Layers] [#GridWidth] [#GridHeight] [#Steps] [#Repetitions]\"." << std::endl;
                std::cout << "Proceeding with default values." << std::endl;
                break;
            case 2 :
                layers = std::stoi(argv[1]);
                break;
            case 3 :
                layers = std::stoi(argv[1]);
                width = std::stoi(argv[2]);
                break;
            case 4 :
                layers = std::stoi(argv[1]);
                width = std::stoi(argv[2]);
                height = std::stoi(argv[3]);
                break;
            case 5 :
                layers = std::stoi(argv[1]);
                width = std::stoi(argv[2]);
                height = std::stoi(argv[3]);
                steps = std::stoi(argv[4]);
                break;
            case 6 :
                layers = std::stoi(argv[1]);
                width = std::stoi(argv[2]);
                height = std::stoi(argv[3]);
                steps = std::stoi(argv[4]);
                reps = std::stoi(argv[5]);
                break;
            default :
                std::cout << "Configuration is available through program arguments in the form \"./WorldEngineTest [#Layers] [#GridWidth] [#Steps] [#Repetitions]\"." << std::endl;
                std::cout << "Discarding parameters beyond " << argv[5] << "." << std::endl;
                layers = std::stoi(argv[1]);
                width = std::stoi(argv[2]);
                height = std::stoi(argv[3]);
                steps = std::stoi(argv[4]);
                reps = std::stoi(argv[5]);
                break;

        }
//	std::cout << "\nTesting class \033[1;34m'worldEngine'\033[0m:" << std::endl;
	
//	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
    testConstructors();
	
//	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
    testGetters();
	
//	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
//    for (int i = 0; i < reps; ++i)
        testMiscellany(layers, width, height, steps, 30000, 30000, 10000);
    /*for (int i = 0; i < 10; ++i)
        testMiscellany(14, 512, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(17, 512, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(20, 512, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(3, 1024, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(5, 1024, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(8, 1024, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(11, 1024, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(14, 1024, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(17, 1024, 500);
    for (int i = 0; i < 10; ++i)
        testMiscellany(20, 1024, 500);*/
//	std::cout << "\nTesting of class \033[1;34m'worldEngine'\033[0m complete." << std::endl;
}
