#include "advectionEngine.h"
#include "airLayer.h"
#include <algorithm>
#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "planet.h"
#include <vector>

typedef double valType;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
using namespace PWM::Engine;

void testConstructors(){
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto l = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 100, 100, 512));
	auto advEn = advectionEngine<dsType, valType>(512, 60.0f, true, true);
	advEn.addLayer(l);
	
	auto x = advEn.getAirLayers();
	if (std::find(x.begin(), x.end(), l) != x.end())
        std::cout << "\033[1;32mLayers end up in the correct vector structure.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Layers not in correct vector structure!\033[0m" << std::endl;
}

void testComparators(){
	
}

void testGetters(){
	auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto l = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 100, 100, 512));
	auto advEn = advectionEngine<dsType, valType>(512, 60.0f, true, true);
	advEn.addLayer(l);
	
	auto x = advEn.getAirLayers();
	if (x.size() == 1 && x[0] == l)
		std::cout << "\033[1;32mgetAirLayers() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAirLayers() does not work as expected!\033[0m" << std::endl;
	
	if (advEn.getDt() == 60.f)
		std::cout << "\033[1;32mgetDt() works as expected.\033[0m" << std::endl;
	else
        std::cout << "\033[1;31mError! getDt() does not work as expected!\033[0m" << std::endl;
}

void print(const dsType& a){
	for (int i = 0; i < a.getWidth(); ++i){
		std::string s = std::to_string(a.getData(i, 0));
		for (int j = 1; j < a.getWidth(); ++j)
			s += "\t" + std::to_string(a.getData(i, j));
		std::cout << s << std::endl;
	}
}

void testSetters(){
	
}

void testMiscellany(int wCoriolisSteps, int noCoriolisSteps, int latCells){
	int threads = omp_get_max_threads();
	auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto l1 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 100, 100, latCells));
	auto advEn = advectionEngine<dsType, valType>(latCells, 1.0f, true, true);
	auto advEnNoCor = advectionEngine<dsType, valType>(latCells, 1.0f, true, false);
	l1->randomInit();
	for (int i = 0; i < l1->getVelocityTheta().size(); ++i){
		l1->setVelocityTheta(i, i);
		l1->setVelocityPhi(i, i);
	}
	auto l2 = l1;
	auto lUnchanged = l1;
	advEn.addLayer(l1);
	advEnNoCor.addLayer(l2);
	
	//std::cout << "Before advection:" << std::endl;
	/*print((*l1).getPressure());
	/*std::cout << std::endl;
	print((*l1).getVelocityTheta());
	std::cout << std::endl;
	print((*l1).getVelocityPhi());*/
	for (int i = 0; i < wCoriolisSteps; ++i){
		advEn.step();
		/*std::cout << "\nAfter step " << i + 1 << ":" << std::endl;
		print((*l1).getPressure());
		std::cout << std::endl;
		print((*l1).getVelocityTheta());
		std::cout << std::endl;
		print((*l1).getVelocityPhi());*/
	}
	std::cout << "\033[1;37mExecution time using " << threads << " threads for " << wCoriolisSteps << " steps of " << latCells << " x " << latCells << " cells with coriolis is: \033[1;33m" << advEn.getRunTimePassed() << " \033[1;37mseconds." << std::endl; 
	if (*l1 != *lUnchanged)
		std::cout << "\033[1;32mAdvection gives some result. Whether those results are correct remains to be seen.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Advection step results in no changes!\033[0m" << std::endl;

	for (int i = 0; i < noCoriolisSteps; ++i)
		advEnNoCor.step();
	std::cout << "\033[1;37mExecution time using " << threads << " threads for " << noCoriolisSteps << " steps of " << latCells << " x " << latCells << " cells without coriolis is: \033[1;33m" << advEnNoCor.getRunTimePassed() << " \033[1;37mseconds." << std::endl; 
	if (*l1 != *l2)
		std::cout << "\033[1;32mCoriolis step results in differences.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Coriolis step makes no difference!\033[0m" << std::endl;
}

int writeTempImage(std::string file, PWM::PWMDataStructure::square2DArray<valType> data){
    std::ofstream fil;
    try {
        fil.open(file);
        //Header
        fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

        //Image data
        auto mM = std::make_pair<double, double>(200, 400);
        double spread = mM.second - mM.first;
        double max = 1200;
        for (int i = 0; i < data.size(); ++i){
            int r, g, b;
            float tNormed, t = data.getData(i);
            if (t < mM.second){
                tNormed = (std::abs(spread) < 1 ? 0 : ((t - mM.first) / spread) * 255);
                r = 255 - (255 - tNormed);
                g = 0;
                b = 255 - tNormed;
            }
            else{
                tNormed = (((t - mM.second) / (max - mM.second)) * 255);
                r = 255;
                g = 255 - (255 - tNormed);
                b = 255 - (255 - tNormed);
            }
            fil << r << " " << g << " " << b << " ";
        }

        fil.close();
        return 0;
    } catch (std::exception& e) {
        std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
        return -1;
    }
}

int writeMoisImage(std::string file, std::shared_ptr<PWM::Model::planet> planetos, PWM::PWMDataStructure::square2DArray<valType> data){
    std::ofstream fil;
    try {
        fil.open(file);
        //Header
        fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

        //Image data
        auto mM = std::make_pair<double, double>(0, 5);
        double spread = mM.second - mM.first;
        for (int i = 0; i < data.size(); ++i){
            float moisNormed = (std::abs(spread) < 1 ? 0 : ((data.getData(i) - mM.first) / spread) * 255);
            int r = 255 - (255 - moisNormed);
            int g = 0;
            int b = 255 - moisNormed;
            fil << r << " " << g << " " << b << " ";
        }

        fil.close();
        return 0;
    } catch (std::exception& e) {
        std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
        return -1;
    }
}

void runImageTest(){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("AdvectionTest");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
            std::filesystem::remove(dir_Entry.path());
    }
    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto l = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 250, 500, 128, 10000.f));
    for (int j = 0; j < l->getPressure().size(); ++j){
        l->setVelocityTheta(j, -0.5);
        l->setVelocityPhi(j, 0.5);
    }
    int startj = 40, startk = 40;
    valType highT = 380, highM = 4;
    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < 10; ++k){
                l->setTemperature(j + startj, k + startk, highT);
                l->setCondensedWater(j + startj, k + startk, highM);
            }

        if ((i + 1) % 2)
            startk += 40;
        else{
            startj += 40;
            startk = 40;
        }
    }

    auto z = advectionEngine<dsType, valType>(128, 10000.0, 60.f, true);
    z.addLayer(l);
    for (int i = 0; i < 300; ++i){
//        if (i % 500 == 0){
                std::stringstream t, c;
                t << "../ImageOutput/AdvectionTest/Temperature_Step_" << i << ".ppm";
                writeTempImage(t.str(), l->getTemperature());
                c << "../ImageOutput/AdvectionTest/Cloud_Step_" << i << ".ppm";
                writeMoisImage(c.str(), planetos, l->getCondensedWater());
//        }
        z.step();
    }
}

int main(int argc, char** argv){
	std::cout << "Testing Eigen inclusion:" << std::endl;
	Eigen::MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 2.5;
	m(0,1) = -1;
	m(1,1) = m(1,0) + m(0,1);
	std::cout << m << std::endl;

	std::cout << "\nTesting class \033[1;34m'AdvectionEngine'\033[0m:" << std::endl;
	int corSteps = 10, noCorSteps = 10, latCells = 512;
	switch (argc){
		case 1 :
			std::cout << "Configuration is available through program arguments in the form \"./AdvectionEngineTest [StepsWithCoriolis] [StepsWithoutCoriolis] [#CellsPerSide]\"." << std::endl;
			std::cout << "Proceeding with default values." << std::endl;
			break;
		case 2 :
			corSteps = std::stoi(argv[1]);
			break;
		case 3 :
			corSteps = std::stoi(argv[1]);
			noCorSteps = std::stoi(argv[2]);
			break;
		case 4 :
			corSteps = std::stoi(argv[1]);
			noCorSteps = std::stoi(argv[2]);
			latCells = std::stoi(argv[3]);
			break;
		default :
			std::cout << "Configuration is available through program arguments in the form \"./AdvectionEngineTest [StepsWithCoriolis] [StepsWithoutCoriolis] [#LatitudeCells]\"." << std::endl;
			std::cout << "Discarding parameters beyond " << argv[3] << "." << std::endl;
			corSteps = std::stoi(argv[1]);
			noCorSteps = std::stoi(argv[2]);
			latCells = std::stoi(argv[3]);
			break;

	}

	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
	testGetters();
	
	std::cout << "\n\033[1;37mTesting advection for " << corSteps << " steps with the Coriolis effect, " << noCorSteps << " steps without the Coriolis effect, both over a grid of " << latCells << " x " << latCells << " cells:\033[0m" << std::endl;
    testMiscellany(corSteps, noCorSteps, latCells);

    std::cout << "\n\033[1;37mTesting advection and writing images to ../ImageOutput/AdvectionTest/. High temp and cloudy regions should move linearly.\033[0m" << std::endl;
    runImageTest();
	
	std::cout << "\nTesting of class \033[1;34m'AdvectionEngine'\033[0m complete." << std::endl;
}
