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
#include "planet.h"
#include "pressureEngine.h"
#include "vizUtils.h"

typedef double valType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType2;

void testConstructors(){

}

void testComparators(){

}

void testGetters(){

}

void testSetters(){

}

void testMiscellany(const size_t nX, const size_t nY, const int steps, const int frac){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/FlatMACUlyssesComparison/FlatMAC");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && ((dir_Entry.path().extension().string() == ".png") || (dir_Entry.path().extension().string() == ".ppm")))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir.remove_filename().concat("Ulysses");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && ((dir_Entry.path().extension().string() == ".png") || (dir_Entry.path().extension().string() == ".ppm")))
            std::filesystem::remove(dir_Entry.path());
    }

    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto l1 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 100, 100, nX, nY, 10000, 2500));
    auto l2 = std::make_shared<PWM::Model::airLayer<dsType2, valType>>(PWM::Model::airLayer<dsType2, valType>(planetos, 100, 100, nX, 10000));
    l1->getObstacles().copy(0);
    l2->getObstacles().copy(0);
    int UlyssesjMax = (nX / 2) - (nY / 2) + 1;
    int UlysseskMin = (nX / 2) + (nY / 2) - 1;
    for (int i = 0; i < nX; ++i){//For width
        l1->setObstacles(0, i, 1);
        l1->setObstacles(nY - 1, i, 1);
        for (int j = 0; j < UlyssesjMax; ++j){
            l2->setObstacles(j, i, 1);
        }
        for (int k = UlysseskMin; k < nX; ++k){
            l2->setObstacles(k, i, 1);
        }
    }
    auto radius = nY / 8;
    std::pair<size_t, size_t> hillCentre = std::make_pair(nY / 2, nX / 5);
    int iMin = hillCentre.first - radius, iMax = hillCentre.first + radius;
    int jMin = hillCentre.second - radius, jMax = hillCentre.second + radius;
    int UlyssesiAdjustment = (nX / 2) - (nY / 2);
    for (int i = iMin; i < iMax; ++i){
        for (int j = jMin; j < jMax; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(hillCentre, std::make_pair(i, j));
            if (dist <= radius){
                l1->setObstacles(i, j, 1.);
                l2->setObstacles(i + UlyssesiAdjustment, j, 1.);
            }
        }
    }
    std::stringstream fO, fO1;
    fO << "../ImageOutput/FlatMACUlyssesComparison/FlatMAC/Layer_Obstacles.ppm";
    fO1 << "../ImageOutput/FlatMACUlyssesComparison/Ulysses/Layer_Obstacles.ppm";
    PWM::Utils::writeTerrElevImage(fO.str(), l1->getObstacles());
    PWM::Utils::writeTerrElevImage(fO1.str(), l2->getObstacles());
//    l1->getVelocityPhi().randomInit(-20, 20);
//    l1->getVelocityTheta().randomInit(-20, 20);
    l1->getVelocityPhi().copy(10);
    l1->getVelocityTheta().copy(0);
    l2->getVelocityPhi().copy(10);
    l2->getVelocityTheta().copy(0);
    int hotWidth = nY / 10;
    int startj = nY / 2, startk = 0;
    double highT = 620.0;
    for (int i = -hotWidth; i < hotWidth; ++i){
        for (int j = -hotWidth; j < hotWidth; ++j){
            l1->setTemperature(i + startj, j + startk, highT);
            l2->setTemperature(i + startj + UlyssesiAdjustment, j + startk, highT);
        }
    }

    auto x = PWM::Engine::advectionEngine<dsType, valType>(nX, nY, 10000, 2500, 1.0, true, false);
    auto xP = PWM::Engine::pressureEngine<dsType, valType>(nX, nY, 10000, 2500, 1.0, true);
    x.addLayer(l1);
    xP.addLayer(l1);
    auto y = PWM::Engine::advectionEngine<dsType2, valType>(nX, 10000, 1.0, true, false);
    auto yP = PWM::Engine::pressureEngine<dsType2, valType>(nX, 10000, 1.0, true);
    y.addLayer(l2);
    yP.addLayer(l2);

    double previousStep = 0;
    /*for (int i = 0; i <= steps; ++i){
        if (i % frac == 0){
            std::stringstream fT, fVY, fVX, fP;
            fT << "../ImageOutput/FlatMACUlyssesComparison/FlatMAC/Layer_Temperature_Step_" << i + 1 << ".ppm";
            fVX << "../ImageOutput/FlatMACUlyssesComparison/FlatMAC/Layer_VelX_Step_" << i + 1 << ".ppm";
            fVY << "../ImageOutput/FlatMACUlyssesComparison/FlatMAC/Layer_VelY_Step_" << i + 1 << ".ppm";
            fP << "../ImageOutput/FlatMACUlyssesComparison/FlatMAC/Layer_Pressure_Step_" << i + 1 << ".ppm";
            PWM::Utils::writeTempImage(fT.str(), l1->getTemperature());
            PWM::Utils::writeVelImage(fVX.str(), l1->getVelocityTheta());
            PWM::Utils::writeVelImage(fVY.str(), l1->getVelocityPhi());
            PWM::Utils::writeTerrElevImage(fP.str(), l1->getPressure());
        }
        previousStep = x.getRunTimePassed();
        x.step();
        xP.step();
        double currentStep = x.getRunTimePassed();
        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;
    }*/

    previousStep = 0;
    for (int i = 0; i <= steps; ++i){
        if (i % frac == 0){
            std::stringstream fT, fVY, fVX, fP;
            fT << "../ImageOutput/FlatMACUlyssesComparison/Ulysses/Layer_Temperature_Step_" << i + 1 << ".ppm";
            fVX << "../ImageOutput/FlatMACUlyssesComparison/Ulysses/Layer_VelX_Step_" << i + 1 << ".ppm";
            fVY << "../ImageOutput/FlatMACUlyssesComparison/Ulysses/Layer_VelY_Step_" << i + 1 << ".ppm";
            fP << "../ImageOutput/FlatMACUlyssesComparison/Ulysses/Layer_Pressure_Step_" << i + 1 << ".ppm";
            PWM::Utils::writeTempImage(fT.str(), l2->getTemperature());
            PWM::Utils::writeVelImage(fVX.str(), l2->getVelocityTheta());
            PWM::Utils::writeVelImage(fVY.str(), l2->getVelocityPhi());
            PWM::Utils::writeTerrElevImage(fP.str(), l2->getPressure());
        }
        previousStep = y.getRunTimePassed();
        y.step();
        yP.step();
        double currentStep = y.getRunTimePassed();
        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;
    }
}

int main(int argc, char** argv){
    std::cout << "\nTesting advection with \033[1;34m'MAC grid'\033[0m:" << std::endl;

    std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
    testMiscellany(1024, 256, 800, 2);

    std::cout << "\nTesting advection with \033[1;34m'MAC grid'\033[0m complete." << std::endl;
}
