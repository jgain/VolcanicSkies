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
