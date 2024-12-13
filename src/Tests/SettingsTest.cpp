#include <fstream>
#include <iostream>
#include "settings.h"

using namespace PWM::Utils;

void testMiscellany(){
    auto p = PWM::Model::planet("../resources/PlanetEarth.json");
    auto s = PWM::Model::sun<float>("../resources/SunSol.json");

    auto x = settings();
    x.P = p;
    x.suns.push_back(s);
    x.terrainFile = "../resources/sthelens_detailed_sub5.obj";
    x.gridWidth = 60;
    
    std::string sFile = "../resources/StandardSettings.json";
    
    auto y = settings(sFile);
    if (x == y)
        std::cout << "\033[1;32mSettings can be read from JSON file.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Settings read from JSON file do not perform as expected!\033[0m" << std::endl;
}

int main(int argc, char** argv){
    std::cout << "\nTesting class \033[1;34m'settings'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'settings'\033[0m complete." << std::endl;
}
