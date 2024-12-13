/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
 */

#include "advectionEngine.h"
#include "airLayer.h"
#include <filesystem>
#include "flatStaggeredGrid.h"
#include <iostream>
#include "obstaclesEngine.h"
#include "pressureEngine.h"
#include "terrain.h"
#include "verticalCouplingEngine.h"
#include "vizUtils.h"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;
using namespace PWM::Engine;

void testGetters(){
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto a = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 1659., 100, 128, 128, 10000, 10000));
    auto t = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, 128, 128, 10000, 10000);
    auto ter = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t);
    valType val = 2315;
    ter->setTemperature(1659, val);

    auto x = obstaclesEngine<dsType, dsSType, valType, valType2>(1., true);
    x.addAirLayer(a);
    x.setTerrain(ter);
    if (x.getAirLayers().at(0)->getHeight() == 1659.)
        std::cout << "\033[1;32mgetAirLayers() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAirLayers() does not perform as expected!\033[0m" << std::endl;

    if (x.getAirLayer(0)->getHeight() == 1659.)
        std::cout << "\033[1;32mgetAirLayer(const size_t index) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAirLayer(const size_t index) does not perform as expected!\033[0m" << std::endl;

    if (x.getTerrain()->getTemperature(1659) == val)
        std::cout << "\033[1;32mgetTerrain() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrain() does not perform as expected!\033[0m" << std::endl;
}

void testSetters(){
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto a = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 1659., 100, 128, 128, 10000, 10000));
    auto t = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, 128, 128, 10000, 10000);
    auto ter = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t);
    valType val = 2315;
    ter->setTemperature(1659, val);

    auto x = obstaclesEngine<dsType, dsSType, valType, valType2>(1., true);
    x.addAirLayer(a);
    x.setTerrain(ter);
    if (x.getAirLayers().at(0)->getHeight() == 1659.)
        std::cout << "\033[1;32maddAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l does not perform as expected!\033[0m" << std::endl;

    if (x.getTerrain()->getTemperature(1659) == val)
        std::cout << "\033[1;32msetTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> t works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> t does not perform as expected!\033[0m" << std::endl;
}

void testSingleHill(const size_t width, const int steps){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/ObstaclesEngineTest/SingleHill");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    valType temp = 273.15;
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));

    auto aL = std::vector<std::shared_ptr<PWM::Model::airLayer<dsType, valType>>>();
    auto cL = std::vector<std::shared_ptr<PWM::Model::convectionLayer<dsType, valType>>>();
    for (int i = 0; i < 6; ++i){
        auto a = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 50 + i * 100, 100, width, width, 10000, 10000));
        a->getTemperature() = temp;
        a->getVelocityTheta() = 0;
        a->getVelocityPhi() = 20;
        aL.push_back(a);
        if (i < 5){
            auto c = std::make_shared<PWM::Model::convectionLayer<dsType, valType>>(PWM::Model::convectionLayer<dsType, valType>(planetos, 100 + i * 100, i * 100, width, width, 10000, 10000));
            a->getTemperature() = temp;
            a->getVelocityTheta() = 0;
            cL.push_back(c);
        }
    }

    for (int i = 9 * width / 20; i < 11 * width / 20; ++i){
        for (int j = width / 20; j < width / 10; ++j){
            for (auto l : aL){
                l->setTemperature(i, j, 640);
            }
        }
    }

    auto convertXYtoI = [](const size_t i, const size_t j, const size_t width)->size_t{
        return (i * width + j);
    };

    auto t = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, width, width, 10000, 10000);
    auto ter = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t);

    valType maxElev = 599, radius = width / 3;
    std::pair<size_t, size_t> layerCentre = std::make_pair(width / 2, width / 2);
    for (int i = layerCentre.first - radius; i < layerCentre.first + radius; ++i){
        for (int j = layerCentre.second - radius; j < layerCentre.second + radius; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
            if (dist <= radius){
                valType normDist = (valType) dist / (valType) radius;
                valType newElev = PWM::Utils::interpolate(maxElev, 0., normDist);
                ter->setElevation(convertXYtoI(i, j, width), newElev);
                for (auto layer : aL){
                    if (newElev >= layer->getLayerTop())
                        layer->setObstacles(i, j, 1.);
                }
            }
        }
    }

    auto x = obstaclesEngine<dsType, dsSType, valType, valType2>(1.f, true);
    auto y = advectionEngine<dsType, valType>(width, width, 10000, 10000, 1.f, true, true);
    auto w = pressureEngine<dsType, valType>(width, width, 10000, 10000, 1.f, true);
    auto z = verticalCouplingEngine<dsType, valType>(width, width, 10000, 10000, 1.f, true);
    for (auto layer : aL){
        x.addAirLayer(layer);
        y.addLayer(layer);
        w.addLayer(layer);
        z.addAirLayer(layer);
    }
    for (auto c : cL){
        z.addConvectionLayer(c);
    }
    x.setTerrain(ter);

    PWM::Utils::writeTerrElevImage("../ImageOutput/ObstaclesEngineTest/SingleHill/terrainModelElev.ppm", x.getTerrain()->getElevation());

    std::cout << "\033[1;36mRunning ObstaclesEngineTest/SingleHill for " << steps << " steps on a grid of " << width << "x" << width << " cells.\033[0m" << std::endl;
    valType previousStep = 0;
    for (int i = 0; i <= steps; ++i){
        previousStep = x.getRunTimePassed() + y.getRunTimePassed() + w.getRunTimePassed() + z.getRunTimePassed();
        if (i % 10 == 0){
            for (int j = 0; j < aL.size(); ++j){
                std::stringstream fT, fsX, fsY;
                fT << "../ImageOutput/ObstaclesEngineTest/SingleHill/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_Temperature_Step_" << i << ".ppm";
                fsX << "../ImageOutput/ObstaclesEngineTest/SingleHill/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_VelocityX_Step_" << i << ".ppm";
                fsY << "../ImageOutput/ObstaclesEngineTest/SingleHill/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_VelocityY_Step_" << i << ".ppm";
                PWM::Utils::writeTempImage(fT.str(), aL[j]->getTemperature());
                PWM::Utils::writeVelImage(fsX.str(), aL[j]->getVelocityPhi());
                PWM::Utils::writeVelImage(fsY.str(), aL[j]->getVelocityTheta());
                if (j < cL.size()){
                    std::stringstream fc;
                    fc << "../ImageOutput/ObstaclesEngineTest/SingleHill/Vertical_Layer_" << j << "_Vertical_Velocity_Step_" << i << ".ppm";
                    PWM::Utils::writeVelImage(fc.str(), cL[j]->getVerticalVelocities());
                }
            }
        }
        y.step();
        x.step();
        w.calcPressures();
//        if (i % 2 == 0){
//            for (int j = 0; j < aL.size(); ++j){
//                std::stringstream fP;
//                fP << "../ImageOutput/ObstaclesEngineTest/SingleHill/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_Pressure_Step_" << i << "_Before_Vert.ppm";
//                PWM::Utils::writeTerrElevImage(fP.str(), aL[j]->getPressure());
//            }
//        }
        z.step();
        w.solvePressureProjection();
        w.step();
//        if (i % 2 == 0){
//            for (int j = 0; j < aL.size(); ++j){
//                std::stringstream fP;
//                fP << "../ImageOutput/ObstaclesEngineTest/SingleHill/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_Pressure_Step_" << i << "_ZAfter_Vert.ppm";
//                PWM::Utils::writeTerrElevImage(fP.str(), aL[j]->getPressure());
//            }
//        }
//        y.smoothPressure();

        double currentStep = x.getRunTimePassed() + y.getRunTimePassed() + w.getRunTimePassed() + z.getRunTimePassed();
        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;

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

    auto totalT = x.getRunTimePassed() + y.getRunTimePassed() + w.getRunTimePassed() + z.getRunTimePassed();
    std::cout << "Running " << steps << " of a grid of " << width << "x" << width << " cells took " << totalT << " seconds." << std::endl;
    std::cout << "Advection took " << y.getRunTimePassed() << " seconds, pressure solve took " << w.getRunTimePassed() << " seconds, obstacles engine took " << x.getRunTimePassed() << " seconds, and vertical coupling took " << z.getRunTimePassed() << " seconds." << std::endl;
}

void testBernoulliHills(const size_t width, const int steps){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/ObstaclesEngineTest/BernoulliPrinciple");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    valType temp = 273.15;
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));

    auto aL = std::vector<std::shared_ptr<PWM::Model::airLayer<dsType, valType>>>();
    auto cL = std::vector<std::shared_ptr<PWM::Model::convectionLayer<dsType, valType>>>();
    for (int i = 0; i < 6; ++i){
        auto a = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 50 + i * 100, 100, width, width, 10000, 10000));
        a->getTemperature() = temp;
        a->getVelocityTheta() = 0;
        a->getVelocityPhi() = 12;
        aL.push_back(a);
        if (i < 5){
            auto c = std::make_shared<PWM::Model::convectionLayer<dsType, valType>>(PWM::Model::convectionLayer<dsType, valType>(planetos, 100 + i * 100, i * 100, width, width, 10000, 10000));
            cL.push_back(c);
        }
    }

    for (int i = 9 * width / 20; i < 11 * width / 20; ++i){
        for (int j = width / 20; j < width / 10; ++j){
            for (auto l : aL){
                l->setTemperature(i, j, 640);
            }
        }
    }

    auto convertXYtoI = [](const size_t i, const size_t j, const size_t width)->size_t{
        return (i * width + j);
    };

    auto t = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, width, width, 10000, 10000);
    auto ter = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t);

    valType maxElev = 599, radius = width / 4;
    std::pair<size_t, size_t> hillCentre = std::make_pair(width / 4, width / 2);
    for (int i = hillCentre.first - radius; i < hillCentre.first + radius; ++i){
        for (int j = hillCentre.second - radius; j < hillCentre.second + radius; ++j){
            valType dist = std::max(std::abs((valType) hillCentre.first - i), std::abs((valType) hillCentre.second - j));
            if (dist <= radius){
                valType normDist = (valType) dist / (valType) radius;
                valType newElev = PWM::Utils::interpolate(maxElev, 0., normDist);
                ter->setElevation(convertXYtoI(i, j, width), newElev);
                for (auto layer : aL){
                    if (newElev >= layer->getLayerTop())
                        layer->setObstacles(i, j, 1.);
                }
            }
        }
    }

    hillCentre = std::make_pair(3 * (width / 4), width / 2);
    for (int i = hillCentre.first - radius; i < hillCentre.first + radius; ++i){
        for (int j = hillCentre.second - radius; j < hillCentre.second + radius; ++j){
            valType dist = std::max(std::abs((valType) hillCentre.first - i), std::abs((valType) hillCentre.second - j));
            if (dist <= radius){
                valType normDist = (valType) dist / (valType) radius;
                valType newElev = PWM::Utils::interpolate(maxElev, 0., normDist);
                ter->setElevation(convertXYtoI(i, j, width), newElev);
                for (auto layer : aL){
                    if (newElev >= layer->getLayerTop())
                        layer->setObstacles(i, j, 1.);
                }
            }
        }
    }

    auto x = obstaclesEngine<dsType, dsSType, valType, valType2>(1.f, true);
    auto y = advectionEngine<dsType, valType>(width, width, 10000, 10000, 1.f, true, true);
    auto w = pressureEngine<dsType, valType>(width, width, 10000, 10000, 1.f, true);
    auto z = verticalCouplingEngine<dsType, valType>(width, width, 10000, 10000, 1.f, true);
    for (auto layer : aL){
        x.addAirLayer(layer);
        y.addLayer(layer);
        w.addLayer(layer);
        z.addAirLayer(layer);
    }
    for (auto c : cL){
        z.addConvectionLayer(c);
    }
    x.setTerrain(ter);

    PWM::Utils::writeTerrElevImage("../ImageOutput/ObstaclesEngineTest/BernoulliPrinciple/terrainModelElev.ppm", x.getTerrain()->getElevation());

    std::cout << "\033[1;36mRunning ObstaclesEngineTest/BernoulliPrinciple for " << steps << " steps on a grid of " << width << "x" << width << " cells.\033[0m" << std::endl;
    valType previousStep = 0;
    for (int i = 0; i <= steps; ++i){
        previousStep = x.getRunTimePassed() + y.getRunTimePassed() + w.getRunTimePassed() + z.getRunTimePassed();
        if (i % 10 == 0){
            for (int j = 0; j < aL.size(); ++j){
                std::stringstream fT, fsX, fsY;
                fT << "../ImageOutput/ObstaclesEngineTest/BernoulliPrinciple/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_Temperature_Step_" << i << ".ppm";
                fsX << "../ImageOutput/ObstaclesEngineTest/BernoulliPrinciple/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_VelocityX_Step_" << i << ".ppm";
                fsY << "../ImageOutput/ObstaclesEngineTest/BernoulliPrinciple/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_VelocityY_Step_" << i << ".ppm";
                PWM::Utils::writeTempImage(fT.str(), aL[j]->getTemperature());
                PWM::Utils::writeVelImage(fsX.str(), aL[j]->getVelocityPhi());
                PWM::Utils::writeVelImage(fsY.str(), aL[j]->getVelocityTheta());
                if (j < cL.size()){
                    std::stringstream fc;
                    fc << "../ImageOutput/ObstaclesEngineTest/BernoulliPrinciple/Vertical_Layer_" << j << "_Vertical_Velocity_Step_" << i << ".ppm";
                    PWM::Utils::writeVelImage(fc.str(), cL[j]->getVerticalVelocities());
                }
            }
        }
        y.step();
        x.step();
        w.calcPressures();
//        if (i % 10 == 0){
//            for (int j = 0; j < aL.size(); ++j){
//                std::stringstream fP;
//                fP << "../ImageOutput/ObstaclesEngineTest/SingleHill/Layer_" << j << "_(" << aL[j]->getHeight() << "m)_Pressure_Step_" << i << ".ppm";
//                PWM::Utils::writeTerrElevImage(fP.str(), aL[j]->getPressure());
//            }
//        }
        z.step();
        w.solvePressureProjection();
        w.step();
//        y.smoothPressure();

        double currentStep = x.getRunTimePassed() + y.getRunTimePassed() + w.getRunTimePassed() + z.getRunTimePassed();
        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;

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

    auto totalT = x.getRunTimePassed() + y.getRunTimePassed() + w.getRunTimePassed() + z.getRunTimePassed();
    std::cout << "Running " << steps << " of a grid of " << width << "x" << width << " cells took " << totalT << " seconds." << std::endl;
    std::cout << "Advection took " << y.getRunTimePassed() << " seconds, pressure solve took " << w.getRunTimePassed() << " seconds, obstacles engine took " << x.getRunTimePassed() << " seconds, and vertical coupling took " << z.getRunTimePassed() << " seconds." << std::endl;
}

int main(int argc, char** argv){
    std::cout << "\nTesting class \033[1;34m'obstaclesEngine'\033[0m:" << std::endl;

    std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
    testGetters();

    std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
    testSetters();

    std::cout << "\n\033[1;37mTesting single hill obstacle:\033[0m" << std::endl;
    testSingleHill(256, 1000);

    std::cout << "\n\033[1;37mTesting double hill for bernoulli principle obstacle:\033[0m" << std::endl;
    testBernoulliHills(256, 1000);

    std::cout << "\nTesting of class \033[1;34m'obstaclesEngine'\033[0m complete." << std::endl;
}
