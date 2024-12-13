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
/*    if (condition)
        std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
*/

#include <iostream>
#include "terrain_structure.hpp"
#include "world.h"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
typedef PWM::PWMDataStructure::square2DArray<valType2> dsSType;
using namespace PWM::Model;

void testConstructors(){
    auto a = world<dsType, dsSType, valType, valType2>();
	bool p1, p2, p3, p4;
	p1 = a.getTerrain() == nullptr;
	p2 = a.getPlanet() == nullptr;
	p3 = a.getAirLayers().size() == 0;
	p4 = a.getConvectionLayers().size() == 0;
	if (p1 && p2 && p3 && p4)
		std::cout << "\033[1;32mDefault constructor works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Default constructor does not perform as expected!\033[0m" << std::endl;

	
	auto b = world<dsType, dsSType, valType, valType2>();
	b.init("../resources/PlanetEarth.json", "../resources/test.pwea");
	p1 = b.getTerrain()->getElevation().size() == 32;
	p2 = b.getPlanet()->getRadius() == 6371000;
	p3 = b.getAirLayers().size() > 0;
	p4 = b.getConvectionLayers().size() == b.getAirLayers().size() - 1;
	if (p1 && p2 && p3 && p4)
		std::cout << "\033[1;32mDefault constructor with init works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Default constructor with init does not perform as expected!\033[0m" << std::endl;

	auto c = world<dsType, dsSType, valType, valType2>("../resources/PlanetEarth.json", "../resources/test.pwea");
	p1 = c.getTerrain()->getElevation().size() == 32;
	p2 = c.getPlanet()->getRadius() == 6371000;
	p3 = c.getAirLayers().size() > 0;
	p4 = c.getConvectionLayers().size() == c.getAirLayers().size() - 1;
	if (p1 && p2 && p3 && p4)
		std::cout << "\033[1;32mDefault constructor with file names works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Default constructor with file names does not perform as expected!\033[0m" << std::endl;
	
	auto s = PWM::Utils::settings("../resources/StandardSettings.json");
	auto d = world<dsType, dsSType, valType, valType2>(s);
	p2 = d.getPlanet()->getRadius() == 6371000;
	p3 = d.getAirLayers().size() > 0;
	p4 = d.getConvectionLayers().size() == d.getAirLayers().size() - 1;
	if (p2 && p3 && p4)
		std::cout << "\033[1;32mConstructor with settings file works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Default constructor with settings file does not perform as expected!\033[0m" << std::endl;
	
	
}

void testComparators(){
	auto a = world<dsType, dsSType, valType, valType2>("../resources/PlanetEarth.json", "../resources/test.pwea");
	auto b = world<dsType, dsSType, valType, valType2>("../resources/PlanetEarth.json", "../resources/test.pwea");
	
	if (a == a)
		std::cout << "\033[1;32mSelf comparison works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Self comparison does not perform as expected!\033[0m" << std::endl;
	
	if (a == b)
		std::cout << "\033[1;32mComparison for identical constructors works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Comparison for identical constructors does not perform as expected!\033[0m" << std::endl;

}

void testGetters(){
	
}

void testSetters(){
	
}

void testMiscellany(valType2 sampleTerFile){
	auto a = world<dsType, dsSType, valType, valType2>("../resources/PlanetEarth.json", "../resources/test.pwea");
	auto p = std::make_shared<planet>("../resources/PlanetEarth.json");
	if (*p == *a.getPlanet())
		std::cout << "\033[1;32minit() gives the expected planet.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! init() does not give the expected planet!\033[0m" << std::endl;
	
	if (a.getTerrain()->getElevation(1) == -15)
		std::cout << "\033[1;32minit() gives the expected terrain elevation.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! init() does not give the expected terrain elevation!\033[0m" << std::endl;
	
	if (a.getTerrain()->getTerrainType(9) == "WATER" && a.getTerrain()->getTerrainType(10) == "SOIL")
		std::cout << "\033[1;32minit() gives the expected terrain type.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! init() does not give the expected terrain type!\033[0m" << std::endl;

	auto b = new world<dsType, dsSType, valType, valType2>("../resources/PlanetEarth.json", sampleTerFile);
	if (b->getTerrain()->getElevation().size() == 1024 * 1024)
		std::cout << "\033[1;32mLoading the Earth map seems to work as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Earth map gives issues!\033[0m" << std::endl;
	delete b;
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'World'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting comparators:\033[0m" << std::endl;
	testComparators();
	
	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
	testGetters();
	
	std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
	testSetters();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany("../terrainMaps/mars/mars.pwea");
	
	std::cout << "\nTesting of class \033[1;34m'World'\033[0m complete." << std::endl;
}
