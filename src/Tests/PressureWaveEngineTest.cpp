#include "cloudUtils.h"
#include "flatStaggeredGrid.h"
#include <filesystem>
#include <iostream>
#include "phaseTransitionEngine.h"
#include "planet.h"
#include "pressureWave.h"
#include "pressureWaveEngine.h"
#include "vdbExporter.h"
#include "vizUtils.h"
#include <vector>

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;
using namespace PWM::Engine;

void testConstructors(){
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto a = PWM::Model::pressureWave<valType>(planetos, 10);

    if (a.getStartTime() == 10)
        std::cout << "\033[1;32mDefault constructor works as expected for start time.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Default constructo does not set start time as expected!\033[0m" << std::endl;

    if (a.getInitialPower() == 1000)
        std::cout << "\033[1;32mDefault constructor works as expected for initial power.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Default constructo does not set initial power as expected!\033[0m" << std::endl;

    if (!a.isActive())
        std::cout << "\033[1;32mDefault constructor works as expected for is active.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Default constructo does not set isActive as expected!\033[0m" << std::endl;
}

void testComparators(){

}

void testGetters(){

}

void testSetters(){

}

void testMiscellany(int width, int steps, valType moisVal){
    auto currDir = std::filesystem::current_path().remove_filename();
    currDir.concat("ImageOutput/PressureWaveEngineTest");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png"))
            std::filesystem::remove(dir_Entry.path());
    }

    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto a = std::make_shared<PWM::Model::pressureWave<valType>>(planetos, 0.8, 100000);
    a->setCentre(5000, 5000);
    a->setHeight(550);

    auto aL = std::vector<std::shared_ptr<PWM::Model::airLayer<dsType, valType>>>();
    for (int i = 0; i < 10; ++i){
        auto l = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 50 + i * 100, 100, width, width, 10000, 10000));
        l->getMoisture() = moisVal;
        l->getCondensedWater() = 0.06;
        l->getTemperature() = 310 - i;
        aL.push_back(l);
    }
    auto cL = std::vector<std::shared_ptr<PWM::Model::convectionLayer<dsType, valType>>>();
    for (int i = 0; i < 9; ++i){
        auto l = std::make_shared<PWM::Model::convectionLayer<dsType, valType>>(PWM::Model::convectionLayer<dsType, valType>(planetos, 50 + i * 100, 150 + i * 100, width, width, 10000, 10000));
        cL.push_back(l);
    }

    auto x = pressureWaveEngine<dsType, valType>(0.5, true);
    x.addExplosion(a);

    auto y = phaseTransitionEngine<dsType, valType>(0.5, true, true, true, true, true, 0.1, 0.1);
    y.setDoWilsonClouds(true);

    valType cloudThreshold = 0.1;
    for (auto l : aL){
        x.addLayer(l);
        y.addAirLayer(l);
    }

    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    double previousStep = 0;
    for (int i = 0; i < steps; ++i){
        if (i == 12)
            std::cerr << "Debug point reached" << std::endl;
        if (i % 1 == 0){
            PWM::Utils::basicCloudExtract(aL, cL, cloudThreshold);
            // vdbExporter.write(aL, cL, "../output", i, x.getSimTimePassed(), cloudThreshold);
            for (int j = 0; j < aL.size(); ++j){
                std::stringstream fC, fM, fCloud;
                fC << "../ImageOutput/PressureWaveEngineTest/Layer_" << j << "_CondensedWater_Step_" << i << ".png";
                fM << "../ImageOutput/PressureWaveEngineTest/Layer_" << j << "_Moisture_Step_" << i << ".png";
                fCloud << "../ImageOutput/PressureWaveEngineTest/Layer_" << j << "_Cloud_Step_" << i << ".png";
                PWM::Utils::writeMoisImage(fC.str(), aL.at(j)->getCondensedWater(), aL.at(j)->getObstacles());
                PWM::Utils::writeMoisImage(fM.str(), aL.at(j)->getMoisture(), aL.at(j)->getObstacles());
                PWM::Utils::writeCloudImage(fCloud.str(), aL.at(j)->getClouds());
            }
        }
        previousStep = x.getRunTimePassed();
        x.step();
        double currentStep = x.getRunTimePassed();
        y.step();
        std::cout << "Step " << i << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
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
}

int main(int argc, char ** argv){
    std::cout << "\nTesting class \033[1;34m'pressureWaveEngine'\033[0m:" << std::endl;

    std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
    testConstructors();

    std::cout << "\n\033[1;37mTesting comparators:\033[0m" << std::endl;
    testComparators();

    std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
    testGetters();

    std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
    testSetters();

    std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
    testMiscellany(512, 100, (argc > 1) ? std::stod(argv[1]) : 0.065);

    std::cout << "\nTesting of class \033[1;34m'pressureWaveEngine'\033[0m complete." << std::endl;
}
