/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/
#include "advectionEngine.h"
#include "airLayer.h"
#include "atmosphericHeatingEngine.h"
#include "convectionLayer.h"
#include <filesystem>
#include "flatStaggeredGrid.h"
#include <fstream>
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
using namespace PWM::Model;

void print(const dsType& a){
	for (int i = 0; i < a.size(); ++i){
		std::string s = "";
		s += "\t" + std::to_string(a.getData(i));
		std::cout << s << std::endl;
	}
}

//Create blank terrain

terrain<dsType, dsSType, valType, valType2> constructTerr(std::shared_ptr<planet> p, int width, float actualSize){
    terrain<dsType, dsSType, valType, valType2> res = terrain<dsType, dsSType, valType, valType2>(p, width, width, actualSize, actualSize);
    valType initVal = 0.f;
    valType2 initStr = "SOIL";
    for (int i = 0; i < res.getElevation().size(); ++i){
        res.setMoisture(i, initVal);
        res.setTerrainType(i, initStr);
    }
    return res;
}

void testConstructors(){
    auto en = verticalCouplingEngine<dsType, valType>(60.f, true, true);
	if (en.getDt() == 60.f)
 		std::cout << "\033[1;32mConstructor runs without dying, i.e., works as expected.\033[0m" << std::endl;
 	else
		std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
}

void testComparators(){
	
}

void testGetters(){
	
}

void testSetters(){
	
}

void testMiscellany(){
	auto planetos = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
    auto en = verticalCouplingEngine<dsType, valType>(512, 20000.f, 60.f, true, true);
    auto l1 = std::make_shared<airLayer<dsType, valType>>(airLayer<dsType, valType>(planetos, 500, 1000, 512, 512, 20000.f, 20000.f));
    auto l2 = std::make_shared<airLayer<dsType, valType>>(airLayer<dsType, valType>(planetos, 1500, 1000, 512, 512, 20000.f, 20000.f));

    auto c1 = std::make_shared<convectionLayer<dsType, valType>>(convectionLayer<dsType, valType>(planetos, 1500, 500, 512, 512, 20000.f, 20000.f));

	en.addAirLayer(l1);
	en.addAirLayer(l2);
	en.addConvectionLayer(c1);

	/*std::cout << "Before convection, convection strength is:" << std::endl;
	print(c1->getVerticalVelocities());

	std::cout << "Before convection, top layer temperature is:" << std::endl;
	print(l2->getTemperature());*/
	for (int i = 0; i < 10; ++i)
	en.step();

    std::cout << "Time taken for ten steps of convection for 512x512 grid: " << en.getRunTimePassed() << " seconds." << std::endl;
	/*std::cout << "After convection, convection strength is:" << std::endl;
	print(c1->getVerticalVelocities());

	std::cout << "After convection, top layer temperature is:" << std::endl;
	print(l2->getTemperature());*/
}

void benardCellTest(){
    auto planetos = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/CVEBénardTest/NoCorrection");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/CVEBénardTest/WithAdvectionAndSpeedCorrection");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    auto ter = constructTerr(planetos, 128, 10000.f);
    auto t = std::make_shared<terrain<dsType, dsSType, valType, valType2>>(ter);

    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
    int startj = 40, startk = 40;
    double highT = 360.0;
    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < 10; ++k)
                t->setTemperature(ijToI(128, j + startj, k + startk), highT);

        if ((i + 1) % 2)
            startk += 40;
        else{
            startj += 40;
            startk = 40;
        }
    }
    std::stringstream terTF;
    terTF << "../ImageOutput/CVEBénardTest/Terrain_Init_Temperature.ppm";
    PWM::Utils::writeTempImage(terTF.str(), t->getTemperature());

    std::vector<std::shared_ptr<airLayer<dsType, valType>>> aLs;
    std::vector<std::shared_ptr<convectionLayer<dsType, valType>>> cLs;
    for (int i = 0; i < 10; ++i){
        auto l = std::make_shared<airLayer<dsType, valType>>(airLayer<dsType, valType>(planetos, 500 * (i + 1), 1000 * (i + 1), 128, 128, 10000.f, 10000.f));
        for (int j = 0; j < l->getPressure().size(); ++j){
            l->setVelocityTheta(j, 0);
            l->setVelocityPhi(j, 0);
        }
        aLs.push_back(l);
    }
//    startj = 80;
//    startk = 80;
//    for (int i = 0; i < 4; ++i){
//        for (int j = 0; j < 20; ++j)
//            for (int k = 0; k < 20; ++k)
//                aLs[0]->setTemperature(j + startj, k + startk, highT);

//        if ((i + 1) % 2)
//            startk += 80;
//        else{
//            startj += 80;
//            startk = 80;
//        }
//    }
    aLs[0]->setMoisture(0.2);

    for (int i = 0; i < 9; ++i){
        auto cL = std::make_shared<convectionLayer<dsType, valType>>(convectionLayer<dsType, valType>(planetos, aLs[i]->getHeight(), aLs[i + 1]->getHeight(), 128, 128, 10000.f, 10000.f));
        cLs.push_back(cL);
    }


    auto x = verticalCouplingEngine<dsType, valType>(128, 128, 10000.f, 10000.0, 30.f, true, true);

    auto y = atmosphericHeatingEngine<dsType, dsSType, valType, valType2>(30.f, true);
    y.setTerrain(t);

    auto z = advectionEngine<dsType, valType>(128, 128, 10000.f, 10000.0, 30.f, true);
    auto w = pressureEngine<dsType, valType>(128, 128, 10000.f, 10000.0, 30.f, true);

    for (auto l : aLs){
        x.addAirLayer(l);
        y.addAirLayer(l);
        z.addLayer(l);
        w.addLayer(l);
    }
    for (auto cl : cLs)
        x.addConvectionLayer(cl);

    for (int i = 0; i < 101; ++i){
        if (i % 10 == 0){
            for (int j = 0; j < aLs.size(); ++j){
                std::stringstream s;
                s << "../ImageOutput/CVEBénardTest/WithAdvectionAndSpeedCorrection/Air_Layer_" << j << "_Step_" << i << "_Temperature.ppm";
                PWM::Utils::writeTempImage(s.str(), aLs[j]->getTemperature());
            }
            for (int j = 0; j < cLs.size(); ++j){
                std::stringstream s;
                s << "../ImageOutput/CVEBénardTest/WithAdvectionAndSpeedCorrection/Convection_Layer_" << j << "_Step_" << i << "_VerticalVelocity.ppm";
                PWM::Utils::writeVelImage(s.str(), cLs[j]->getVerticalVelocities());
            }
        }
        y.step();
        z.step();
        w.calcPressures();
        x.step();
        w.solvePressureProjection();
        w.step();
    }
}

void minimalistBenardCellTest(int layers, int steps){
    auto planetos = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/CVEBénardTest");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/CVEBénardTest/MinimalTests");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }

    auto ter = constructTerr(planetos, 512, 10000.f);
    auto t = std::make_shared<terrain<dsType, dsSType, valType, valType2>>(ter);

    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
    int startj = 256, startk = 256;
    int hotWidth = 512 / 15;
    double highT = 660.0;
//    for (int j = -hotWidth; j < hotWidth; ++j)
//        for (int k = -hotWidth; k < hotWidth; ++k)
//            t->setTemperature(ijToI(512, j + startj, k + startk), highT);

    std::pair<size_t, size_t> layerCentre = std::make_pair(256, 256);
    for (int i = layerCentre.first - hotWidth; i < layerCentre.first + hotWidth; ++i){
        for (int j = layerCentre.second - hotWidth; j < layerCentre.second + hotWidth; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
            if (dist <= hotWidth){
                t->setTemperature(ijToI(512, i, j), highT);
            }
        }
    }

    layerCentre = std::make_pair(128, 128);
    for (int i = layerCentre.first - hotWidth; i < layerCentre.first + hotWidth; ++i){
        for (int j = layerCentre.second - hotWidth; j < layerCentre.second + hotWidth; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
            if (dist <= hotWidth){
                t->setTemperature(ijToI(512, i, j), highT);
            }
        }
    }

    layerCentre = std::make_pair(384, 128);
    for (int i = layerCentre.first - hotWidth; i < layerCentre.first + hotWidth; ++i){
        for (int j = layerCentre.second - hotWidth; j < layerCentre.second + hotWidth; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
            if (dist <= hotWidth){
                t->setTemperature(ijToI(512, i, j), highT);
            }
        }
    }

    layerCentre = std::make_pair(128, 384);
    for (int i = layerCentre.first - hotWidth; i < layerCentre.first + hotWidth; ++i){
        for (int j = layerCentre.second - hotWidth; j < layerCentre.second + hotWidth; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
            if (dist <= hotWidth){
                t->setTemperature(ijToI(512, i, j), highT);
            }
        }
    }

    layerCentre = std::make_pair(384, 384);
    for (int i = layerCentre.first - hotWidth; i < layerCentre.first + hotWidth; ++i){
        for (int j = layerCentre.second - hotWidth; j < layerCentre.second + hotWidth; ++j){
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
            if (dist <= hotWidth){
                t->setTemperature(ijToI(512, i, j), highT);
            }
        }
    }

    std::stringstream terTF;
    terTF << "../ImageOutput/CVEBénardTest/Terrain_Init_Temperature.ppm";
    PWM::Utils::writeTempImage(terTF.str(), t->getTemperature());

    std::vector<std::shared_ptr<airLayer<dsType, valType>>> aLs;
    std::vector<std::shared_ptr<convectionLayer<dsType, valType>>> cLs;
    for (int i = 0; i < layers; ++i){
        auto l = std::make_shared<airLayer<dsType, valType>>(airLayer<dsType, valType>(planetos, 250 + 500 * i, 500., 512, 512, 5000.f, 5000.f));
        l->getVelocityTheta() = 0;
        l->getVelocityPhi() = 0;
        l->getMoisture() = 0.1;
        aLs.push_back(l);
    }
//    startj = 80;
//    startk = 80;
//    for (int i = 0; i < 4; ++i){
//        for (int j = 0; j < 20; ++j)
//            for (int k = 0; k < 20; ++k)
//                aLs[0]->setTemperature(j + startj, k + startk, highT);

//        if ((i + 1) % 2)
//            startk += 80;
//        else{
//            startj += 80;
//            startk = 80;
//        }
//    }
    aLs[0]->setMoisture(0.2);

    for (int i = 0; i < layers - 1; ++i){
        auto cL = std::make_shared<convectionLayer<dsType, valType>>(convectionLayer<dsType, valType>(planetos, aLs[i]->getHeight(), aLs[i + 1]->getHeight(), 512, 512, 5000.f, 5000.f));
        cLs.push_back(cL);
    }

    float baseStep = 1.0f;
    float vertEngineStep = 3 * baseStep;

    auto x = verticalCouplingEngine<dsType, valType>(512, 512, 5000.f, 5000.f, vertEngineStep, true, true);

    auto y = atmosphericHeatingEngine<dsType, dsSType, valType, valType2>(baseStep, true);
    y.setTerrain(t);

    auto z = advectionEngine<dsType, valType>(512, 512, 5000.f, 5000.f, baseStep, true, false);
    auto zP = pressureEngine<dsType, valType>(512, 512, 5000.f, 5000.f, baseStep, true);

    auto w = obstaclesEngine<dsType, dsSType, valType, valType2>(baseStep, true);
    w.setTerrain(t);

    for (auto l : aLs){
        w.addAirLayer(l);
        x.addAirLayer(l);
        y.addAirLayer(l);
        zP.addLayer(l);
        z.addLayer(l);
    }
    for (auto cl : cLs)
        x.addConvectionLayer(cl);

    double previousStep = 0;
    double simTime = 0;
    for (int i = 0; i <= steps; ++i){
        if (i % 16 == 0/*> 677 && i < 685*/){
            for (int j = 0; j < aLs.size(); ++j){
                std::stringstream s, s1, s2, s3;
                s << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_Temperature_Before_Step_" << i << ".ppm";
                s1 << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_VelX_Before_Step_" << i << ".ppm";
                s2 << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_VelY_Before_Step_" << i << ".ppm";
//                s3 << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_Pressure_Before_Step_" << i << ".ppm";
                PWM::Utils::writeTempImage(s.str(), aLs[j]->getTemperature());
                PWM::Utils::writeVelImage(s1.str(), aLs[j]->getVelocityTheta());
                PWM::Utils::writeVelImage(s2.str(), aLs[j]->getVelocityPhi());
//                PWM::Utils::writeTerrElevImage(s3.str(), aLs[j]->getPressure());
            }
            for (int j = 0; j < cLs.size(); ++j){
                std::stringstream s;
                s << "../ImageOutput/CVEBénardTest/MinimalTests/Convection_Layer_" << j << "_Before_Step_" << i << "_VerticalVelocity.ppm";
                    PWM::Utils::writeVelImage(s.str(), cLs[j]->getVerticalVelocities());
            }
        }
        previousStep = x.getRunTimePassed() + y.getRunTimePassed() + z.getRunTimePassed() + w.getRunTimePassed();

        y.step();
        z.step();
        w.step();
        zP.step();
        if (!((int) simTime % (int) vertEngineStep)){
            x.step();
            /*//        for (int j = 0; j < aLs.size(); ++j){
//            std::stringstream s, s1, s2, s3;
//            s << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_Temperature_Before_Step_" << i << "_PresSolve.ppm";
//            s1 << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_VelX_Before_Step_" << i << "_PresSolve.ppm";
//            s2 << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_VelY_Before_Step_" << i << "_PresSolve.ppm";
//            s3 << "../ImageOutput/CVEBénardTest/MinimalTests/Air_Layer_" << j << "_" << aLs[j]->getHeight() << "m_Pressure_Before_Step_" << i << "_PresSolve.ppm";
//            PWM::Utils::writeTempImage(s.str(), aLs[j]->getTemperature());
//            PWM::Utils::writeVelImage(s1.str(), aLs[j]->getVelocityTheta());
//            PWM::Utils::writeVelImage(s2.str(), aLs[j]->getVelocityPhi());
//            PWM::Utils::writeTerrElevImage(s3.str(), aLs[j]->getPressure());
//        }
//        for (int j = 0; j < cLs.size(); ++j){
//            std::stringstream s;
//            s << "../ImageOutput/CVEBénardTest/MinimalTests/Convection_Layer_" << j << "_VerticalVelocity_Before_Step_" << i << "_PresSolve.ppm";
//            PWM::Utils::writeVelImage(s.str(), cLs[j]->getVerticalVelocities());
//        }*/
//            zP.solvePressureProjection(); //NB! Not needed with velocity adjust?
        }
//        zP.step();
//        zP.smoothPressure();

        double currentStep = x.getRunTimePassed() + y.getRunTimePassed() + z.getRunTimePassed() + w.getRunTimePassed();
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

        simTime += baseStep;
    }
}

int main(int argc, char** argv){
    std::cout << "\nTesting class \033[1;34m'verticalCouplingEngine'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
//    std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
//	testMiscellany();

//    std::cout << "\n\033[1;37mPerforming Bénard cell test case:\033[0m" << std::endl;
//    benardCellTest();

    std::cout << "\n\033[1;37mPerforming minimal Bénard cell test case (400 steps) with 8 airlayers:\033[0m" << std::endl;
    minimalistBenardCellTest(8, 1600);

    std::cout << "\nTesting of class \033[1;34m'verticalCouplingEngine'\033[0m complete." << std::endl;
}
