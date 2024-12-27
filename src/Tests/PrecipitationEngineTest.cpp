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
/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
 */
#include <iostream>
#include "PerlinNoise.hpp"
#include "square2DArray.h"
#include "vdbExporter.h"
#include "vizUtils.h"
#include "worldEngine.h"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;
using namespace PWM::Engine;

void testConstructors(){

}

void testComparators(){

}

void testGetters(){

}

void testSetters(){

}

void testMiscellany(int width, int steps){
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

    currDir.concat("ImageOutput/RainfallEngineTest");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));

    auto aL = std::vector<std::shared_ptr<PWM::Model::airLayer<dsType, valType>>>();
    auto cL = std::vector<std::shared_ptr<PWM::Model::convectionLayer<dsType, valType>>>();
    for (int i = 0; i < 12; ++i){
        auto a = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 200 + i * 400, 400, width, 10000));
        a->randomInit();
        a->getMoisture() = 0.01;
        a->getCondensedWater() = 0;
        aL.push_back(a);

        if (i < 11){
            auto c = std::make_shared<PWM::Model::convectionLayer<dsType, valType>>(PWM::Model::convectionLayer<dsType, valType>(planetos, 400 + i * 400, i * 400, width, 10000));
            cL.push_back(c);
        }
    }
    aL[3]->getCondensedWater() = 2;
    aL[3]->getMoisture() = 1;

    auto t = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, width, 10000);
    auto ter = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t);
    valType d = 10;
    for (auto j = 0; j < t.getElevation().size(); ++j)
        ter->setElevation(j, d);

    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
    int hotWidth = width / 10;
    int startj = width / 2, startk = width / 2;
    double highT = 800.0;
    for (int j = -hotWidth; j < hotWidth; ++j){
        for (int k = -hotWidth; k < hotWidth; ++k){
            ter->setTemperature(ijToI(width, j + startj, k + startk), highT);
        }
    }

    auto x = precipitationEngine<dsType, dsSType, valType, valType2>(1.0, true, 0.05);
    auto y = advectionEngine<dsType, valType>(width, width, 10000, 10000, 1.0, true, true);
    auto yp = pressureEngine<dsType, valType>(width, width, 10000, 10000, 1.0, true);
    auto z = verticalCouplingEngine<dsType, valType>(width, width, 10000, 10000, 1.0, true, true);

    for (auto i : aL){
        x.addAirLayer(i);
        y.addLayer(i);
        yp.addLayer(i);
        z.addAirLayer(i);
    }

    for (auto i : cL){
        x.addConvectionLayer(i);
        z.addConvectionLayer(i);
    }

    x.setTerrain(ter);

    for (int i = 0; i < steps; ++i){
        if (i % 20 == 0){
            std::stringstream fGW;
            fGW << "../ImageOutput/PrecipitationEngineTest/Rainfall/GroundWater_Step_" << i << ".ppm";
            PWM::Utils::writeMoisImage(fGW.str(), ter->getMoisture());
            for (int j = 0; j < aL.size(); ++j){
                std::stringstream fC, fsX, fsY, fT;
                fC << "../ImageOutput/PrecipitationEngineTest/Rainfall/Air_Layer_" << j << "_(" << aL[j]->getHeight() << "m)_CloudWater_Step_" << i << ".ppm";
                fT << "../ImageOutput/PrecipitationEngineTest/Rainfall/Air_Layer_" << j << "_(" << aL[j]->getHeight() << "m)_Temperature_Step_" << i << ".ppm";
                fsX << "../ImageOutput/PrecipitationEngineTest/Rainfall/Air_Layer_" << j << "_(" << aL[j]->getHeight() << "m)_VelocityX_Step_" << i << ".ppm";
                fsY << "../ImageOutput/PrecipitationEngineTest/Rainfall/Air_Layer_" << j << "_(" << aL[j]->getHeight() << "m)_VelocityY_Step_" << i << ".ppm";
                PWM::Utils::writeMoisImage<valType>(fC.str(), aL[j]->getCondensedWater());
                PWM::Utils::writeTempImage<valType>(fT.str(), aL[j]->getTemperature());
                PWM::Utils::writeVelImage<valType>(fsX.str(), aL[j]->getVelocityPhi());
                PWM::Utils::writeVelImage<valType>(fsY.str(), aL[j]->getVelocityTheta());

                if (j < cL.size()){
                    std::stringstream fc, fR;
                    fc << "../ImageOutput/PrecipitationEngineTest/Rainfall/Vertical_Layer_" << j << "_Vertical_Velocity_Step_" << i << ".ppm";
                    fR << "../ImageOutput/PrecipitationEngineTest/Rainfall/Vertical_Layer_" << j << "_Rainfall_Step_" << i << ".ppm";
                    PWM::Utils::writeVelImage<valType>(fc.str(), cL[j]->getVerticalVelocities());
                    PWM::Utils::writeMoisImage<valType>(fR.str(), cL[j]->getRainfall());
                }
            }
        }
        x.step();
        y.step();
        yp.step();
        z.step();
    }
}

int testAshRain(int steps, int width, int height, valType xSize, valType ySize, int layers, valType topOfAtmo){
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

    currDir.concat("ImageOutput/PrecipitationEngineTest/AshFall");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("PrecipitationEngineTest/output/");
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
    for (int i = 0; i < alSize - 1; ++i){
        y->getAirLayer(i)->getVelocityTheta() = 15 - (0.0025 * (alSize - i));
        y->getAirLayer(i)->getVelocityPhi() = 0 + (0.0025 * (alSize - i));
        y->getAirLayer(i)->getCondensedWater() = 0.01;
        y->getAirLayer(i)->getMoisture() = 0.01;
        y->getAirLayer(i)->getParticulates() = 0;
    }
    y->getAirLayer(layers - 1)->getVelocityTheta() = 0;
    y->getAirLayer(layers - 1)->getVelocityPhi() = 0;
    y->getAirLayer(layers - 1)->getCondensedWater() = 0.01;
    y->getAirLayer(layers - 1)->getMoisture() = 0.01;
    y->getAirLayer(layers - 1)->getParticulates() = 0;

    int ashCloudWidth = width / 4;
    int ashCloudHeight = height / 2;
    valType ashCloudMaxDensity = 4;
    std::pair<size_t, size_t> ashCloudCentre = std::make_pair(height / 2, width / 4);

    const siv::PerlinNoise::seed_type seed = 123456u;
    const siv::PerlinNoise perNoiGen{ seed };

    for (int i = ashCloudCentre.first - ashCloudHeight; i < ashCloudCentre.first + ashCloudHeight; ++i){
        for (int j = ashCloudCentre.second - ashCloudWidth; j < ashCloudCentre.second + ashCloudWidth; ++j){
            bool valid = PWM::Utils::calcPointInEllipse<int, int, int>(ashCloudCentre, std::make_pair(i, j), ashCloudHeight, ashCloudWidth);
            if (valid){
                valType i1 = i + 0.5, j1 = j + 0.5;
                valType ashVal = perNoiGen.normalizedOctave2D_01(i1 * 0.01, j1 * 0.01, 8) * ashCloudMaxDensity;
                y->getAirLayer(layers - 1)->setParticulates(i, j, ashVal);
            }
        }
    }

    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    float weatherDT = 1.f;
    float cveDT = 3.f;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    b.outputDirectory = "../ImageOutput/MergeTest/AshRain";
    auto aL = b.getWorldModel()->getAirLayers();
    auto cL = b.getWorldModel()->getConvectionLayers();

    double previousStep = 0;
    for (int i = 0; i <= steps; ++i){
//        PWM::Utils::writeAirImages("../ImageOutput/PrecipitationEngineTest/AshFall", i, b.getWorldModel());
        if (i % 4 == 0){
            for (int j = 0; j < aL.size(); ++j){
                std::stringstream fA, fPar;
                fPar << "../ImageOutput/PrecipitationEngineTest/AshFall/Air_Layer_" << j << "_" << aL[j]->getHeight() << "m_Particulates_Before_Step_" << i << ".ppm";
                PWM::Utils::writeAshImage(fPar.str(), aL[j]->getParticulates());
                if (j < cL.size()){
                    fA << "../ImageOutput/PrecipitationEngineTest/AshFall/Convection_Layer_" << j << "_Ashfall_Before_Step_" << i << ".ppm";
                    PWM::Utils::writeAshImage(fA.str(), cL[j]->getAshfall());
                }
            }
//            vdbExporter.write(b.getWorldModel(), "../PrecipitationEngineTest/output", i, b.getSimTimePassed());
        }
        previousStep = b.getRunTimePassed();
        b.step();
        double currentStep = b.getRunTimePassed();
        std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
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
    b.printTiming(true);
    return 0;
}

int main(int argc, char** argv){
    std::cout << "\nTesting class \033[1;34m'precipitationEngine'\033[0m:" << std::endl;

//    std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
//    testConstructors();

//    std::cout << "\n\033[1;37mTesting comparators:\033[0m" << std::endl;
//    testComparators();

//    std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
//    testGetters();

//    std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
//    testSetters();

//    std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
//    testMiscellany(128, 500);
    std::cout << "\n\033[1;37mTesting generic ash fall without a plume:\033[0m" << std::endl;
    testAshRain(1000, 512, 512, 15000, 15000, 15, 15000);

    std::cout << "\nTesting of class \033[1;34m'precipitationEngine'\033[0m complete." << std::endl;
}
