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
#include "airLayer.h"
#include "cloudUtils.h"
#include <filesystem>
#include "planet.h"
#include "vdbExporter.h"
#include "world.h"
#include "worldEngine.h"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;

void testMiscellany(int layers, int width){
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

    currDir.remove_filename().remove_filename().concat("ImageOutput/CloudClassificationTests");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir.concat("/PureStratus");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir.remove_filename().concat("SpottedCumulus");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir.remove_filename().concat("StreakedCirrus");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    std::cout << "Creating uniform stratus environment" << std::endl;
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
//    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj");
    auto ter = terrain_structure();
    ter.testTerrainStructure(width);
    y->maxLayers = layers;
    y->init("../resources/PlanetEarth.json", ter);
    y->addSun("../resources/SunSol.json");
    for (int i = 0; i < y->maxLayers; ++i){
        y->getAirLayer(i)->setHeight(100 + i * 2000);
        for (int j = 0; j < y->getAirLayer(i)->getCondensedWater().size(); ++j){
             y->getAirLayer(i)->setCondensedWater(j, 0.2);
        }
    }
    PWM::Utils::basicCloudExtract(y, "../ImageOutput/CloudClassificationTests/PureStratus", 0, 0.1);
    for (int j = 0; j < layers; ++j){
        PWM::Utils::writeCloudImages<valType>("../ImageOutput/CloudClassificationTests/PureStratus", 0, j, y->getAirLayer(j)->getHeight(), y->getAirLayer(j)->getClouds());
    }
    std::cout << "Clouds extracted from uniform stratus environment. See ImageOutput/CloudClassificationTests/PureStratus for output." << std::endl;

    std::cout << "Creating spotted cumulus environment" << std::endl;
    auto z = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    z->maxLayers = layers;
    z->init("../resources/PlanetEarth.json", ter);
    z->addSun("../resources/SunSol.json");
    for (int i = 0; i < layers; ++i){
        z->getAirLayer(i)->setHeight(100 + i * 1500);
        z->getAirLayer(i)->getCondensedWater() = 0;
        if (i == 5){
            int startj = 80, startk = 80;
            for (int m = 0; m < 4; ++m){
                for (int j = 0; j < 20; ++j)
                    for (int k = 0; k < 20; ++k){
                        z->getAirLayer(i)->setCondensedWater(j + startj, k + startk, 8.2);
                        if (i < layers - 1)
                            z->getConvectionLayer(i)->setVerticalVelocity(j + startj, k + startk, 20);
                    }

                if ((m + 1) % 2)
                    startk += 80;
                else{
                    startj += 80;
                    startk = 80;
                }
            }
        }
    }
    PWM::Utils::basicCloudExtract(z, "../ImageOutput/CloudClassificationTests/SpottedCumulus", 0, 0.1);
    for (int j = 0; j < z->maxLayers; ++j){
        PWM::Utils::writeCloudImages<valType>("../ImageOutput/CloudClassificationTests/SpottedCumulus", 0, j, z->getAirLayer(j)->getHeight(), z->getAirLayer(j)->getClouds());
    }
    std::cout << "Clouds extracted from spotted cumulus environment. See ImageOutput/CloudClassificationTests/SpottedCumulus for output." << std::endl;

    std::cout << "Creating streaked cirrus environment" << std::endl;
    auto x = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    x->maxLayers = 5;
    x->init("../resources/PlanetEarth.json", ter);
    x->addSun("../resources/SunSol.json");
    for (int i = 0; i < 5; ++i){
        x->getAirLayer(i)->setHeight(1000 + i * 2000);
        x->getAirLayer(i)->getCondensedWater() = 0;
    }
    int j = x->getAirLayers().size() - 1;
    valType cirrusCloudVal = 0.5, cirrusUpliftVal = 0.325;
    for (int i = 0; i < x->getAirLayer(j)->getObstacles().size(); ++i){
        if (i % 2){
            x->getAirLayer(j)->setCondensedWater(i, cirrusCloudVal);
            x->getConvectionLayer(j - 1)->setVerticalVelocity(i, cirrusUpliftVal);
        }
    }
    PWM::Utils::basicCloudExtract(x, "../ImageOutput/CloudClassificationTests/StreakedCirrus", 0, 0.1);
    for (int j = 0; j < x->maxLayers; ++j){
        PWM::Utils::writeCloudImages<valType>("../ImageOutput/CloudClassificationTests/StreakedCirrus", 0, j, x->getAirLayer(j)->getHeight(), x->getAirLayer(j)->getClouds());
    }
    std::cout << "Clouds extracted from streaked cirrus environment. See ImageOutput/CloudClassificationTests/StreakedCirrus for output." << std::endl;
}

int createTestScene(int width, int height, valType xSize, valType ySize, int layers, valType topOfAtmo){
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

    currDir.concat("ImageOutput/CloudClassificationTests/TestScene");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("CloudClassificationTests/output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }

    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));

    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = terrain_structure();
    ter.testTerrainStructure(width);
    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, width, height, xSize, ySize, layerHeights);

    int alSize = y->getAirLayers().size();
    for (int k = 0; k < alSize; ++k){
        y->getAirLayer(k)->getVelocityTheta() = 50 - (0.005 * (alSize - k));
        y->getAirLayer(k)->getVelocityPhi() = 0;
        y->getAirLayer(k)->getCondensedWater() = 0.03;
        y->getAirLayer(k)->getMoisture() = 0.03;
        y->getAirLayer(k)->getParticulates() = 0;
    }

    int cloudWidth = width / 10;
    int cloudHeight = height / 10;
    valType dense = 0.3;
    valType vertVel = 5;

    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    int widthPortion = width / 4, heightPortion = height / 4;

    //create cumulus

    std::pair<size_t, size_t> cumulusCloudCentre = std::make_pair(heightPortion, widthPortion);
    std::cout << "Cumulonimbus should be in range [" << cumulusCloudCentre.first - cloudHeight << ", " << cumulusCloudCentre.second - cloudWidth << "] to (" << cumulusCloudCentre.first + cloudHeight << ", " << cumulusCloudCentre.second + cloudWidth << "), with the middle at (" << cumulusCloudCentre.first << ", " << cumulusCloudCentre.second << ")." << std::endl;
    valType cumuloTemp = y->getAirLayer(1)->getTemperature(1) + 15;
    for (int i = cumulusCloudCentre.first - cloudHeight; i < cumulusCloudCentre.first + cloudHeight; ++i){
        for (int j = cumulusCloudCentre.second - cloudWidth; j < cumulusCloudCentre.second + cloudWidth; ++j){
            y->getAirLayer(1)->setCondensedWater(i, j, dense);
            y->getAirLayer(1)->setMoisture(i, j, dense * 0.4);
            y->getAirLayer(0)->setTemperature(i, j, cumuloTemp + 10);
            y->getAirLayer(1)->setTemperature(i, j, cumuloTemp);
            y->getConvectionLayer(0)->setVerticalVelocity(i, j, vertVel);
            y->getConvectionLayer(1)->setVerticalVelocity(i, j, vertVel);

        }
    }

    //create stratus
    std::pair<size_t, size_t> stratusCloudCentre = std::make_pair(heightPortion * 3, widthPortion);
    std::cout << "Cumulonimbus should be in range [" << stratusCloudCentre.first - cloudHeight << ", " << stratusCloudCentre.second - cloudWidth << "] to (" << stratusCloudCentre.first + cloudHeight << ", " << stratusCloudCentre.second + cloudWidth << "), with the middle at (" << stratusCloudCentre.first << ", " << stratusCloudCentre.second << ")." << std::endl;
    for (int i = stratusCloudCentre.first - cloudHeight; i < stratusCloudCentre.first + cloudHeight; ++i){
        for (int j = stratusCloudCentre.second - cloudWidth; j < stratusCloudCentre.second + cloudWidth; ++j){
            y->getAirLayer(1)->setCondensedWater(i, j, 0.11);
            y->getAirLayer(1)->setMoisture(i, j, dense * 0.25);
        }
    }

    //create cirrus
    std::pair<size_t, size_t> cirrusCloudCentre = std::make_pair(heightPortion, widthPortion * 3);
    std::cout << "Cumulonimbus should be in range [" << cirrusCloudCentre.first - cloudHeight << ", " << cirrusCloudCentre.second - cloudWidth << "] to (" << cirrusCloudCentre.first + cloudHeight << ", " << cirrusCloudCentre.second + cloudWidth << "), with the middle at (" << cirrusCloudCentre.first << ", " << cirrusCloudCentre.second << ")." << std::endl;
    for (int i = cirrusCloudCentre.first - cloudHeight; i < cirrusCloudCentre.first + cloudHeight; ++i){
        for (int j = cirrusCloudCentre.second - cloudWidth; j < cirrusCloudCentre.second + cloudWidth; ++j){
            y->getAirLayer(layers - 1)->setCondensedWater(i, j, dense * 0.5);
            y->getAirLayer(layers - 1)->setMoisture(i, j, dense * 0.4);
        }
    }

    //create cumulonimbus
    std::pair<size_t, size_t> cumulonimboCloudCentre = std::make_pair(heightPortion * 3, widthPortion * 3);
    std::cout << "Cumulonimbus should be in range [" << cumulonimboCloudCentre.first - cloudHeight << ", " << cumulonimboCloudCentre.second - cloudWidth << "] to (" << cumulonimboCloudCentre.first + cloudHeight << ", " << cumulonimboCloudCentre.second + cloudWidth << "), with the middle at (" << cumulonimboCloudCentre.first << ", " << cumulonimboCloudCentre.second << ")." << std::endl;
    for (int i = cumulonimboCloudCentre.first - cloudHeight; i < cumulonimboCloudCentre.first + cloudHeight; ++i){
        for (int j = cumulonimboCloudCentre.second - cloudWidth; j < cumulonimboCloudCentre.second + cloudWidth; ++j){
            for (int k = 2; k < 9; ++k){
                cumuloTemp = y->getAirLayer(k)->getTemperature(1) + 15;
                y->getAirLayer(k)->setCondensedWater(i, j, dense);
                y->getAirLayer(k)->setMoisture(i, j, dense * 0.4);
                y->getAirLayer(k)->setTemperature(i, j, cumuloTemp);
                y->getConvectionLayer(k - 1)->setVerticalVelocity(i, j, vertVel);
            }
        }
    }

    float weatherDT = 1.f;
    float cveDT = 3.f;
    double previousStep = 0;
    auto b = PWM::Engine::advectionEngine<dsType, valType>(width, height, xSize, ySize, weatherDT, true, false);
    auto bp = PWM::Engine::pressureEngine<dsType, valType>(width, height, xSize, ySize, weatherDT, true);
    auto x = PWM::Engine::verticalCouplingEngine<dsType, valType>(width, height, xSize, ySize, cveDT, true, true);
    for (auto l : y->getAirLayers()){
        b.addLayer(l);
        bp.addLayer(l);
        x.addAirLayer(l);
    }
    for (auto l : y->getConvectionLayers())
        x.addConvectionLayer(l);

    for (int i = 0; i <= 50; ++i){
        PWM::Utils::basicCloudExtract(y, "../ImageOutput/CloudClassificationTests/TestScene", i, 0.1);
        vdbExporter.write(y, "../CloudClassificationTests/output", i, b.getSimTimePassed());
        PWM::Utils::writeAirImages("../ImageOutput/CloudClassificationTests/TestScene", i, y);
        previousStep = b.getRunTimePassed();
        b.step();
        bp.step();
        if (! std::fmod(b.getSimTimePassed(), cveDT))
            x.step();
        double currentStep = b.getRunTimePassed();
        std::cout << "Step " << i + 1 << " of the advection engine took " << currentStep - previousStep << " seconds." << std::endl;
    }
    PWM::Utils::basicCloudExtract(y, "../ImageOutput/CloudClassificationTests/TestScene", 51, 0.1);
    vdbExporter.write(y, "../CloudClassificationTests/output", 51, b.getSimTimePassed());
    PWM::Utils::writeAirImages("../ImageOutput/CloudClassificationTests/TestScene", 51, y);
    return 0;
}

int main(int argc, char** argv){
//	std::cout << "\nTesting \033[1;34m'cloudUtils'\033[0m utilities:" << std::endl;

//	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
//	testConstructors();

//	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
//	testGetters();

    std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
//    for (int i = 0; i < 10; ++i){
//        testMiscellany(7, 256);
//    }

    std::cout << "\n\033[1;37mCreating test scence:\033[0m" << std::endl;
    createTestScene(512, 512, 30000, 30000, 10, 11000);

//	std::cout << "\nTesting of \033[1;34m'cloudUtils'\033[0m utilities complete." << std::endl;
}
