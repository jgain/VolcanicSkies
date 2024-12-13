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
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/

#include "airLayer.h"
#include "heatDissipationEngine.h"
#include <iostream>
#include "planet.h"
#include "square2DArray.h"

typedef double valType;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
using namespace PWM::Engine;

void testConstructors(){
	auto a = heatDissipationEngine<dsType, valType>(60.f, true, 0.001);
	if (a.getCoefficient() != 0)
		std::cout << "\033[1;32mConstructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Constructor does not perform as expected!\033[0m" << std::endl;
}

void testGetters(){
	auto a = heatDissipationEngine<dsType, valType>(60.f, true, 0.02);
	if (a.getCoefficient() == 0.02)
		std::cout << "\033[1;32mgetCoefficient() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getCoefficient() does not perform as expected!\033[0m" << std::endl;
}

void testSetters(){
	auto a = heatDissipationEngine<dsType, valType>(60.f, true, 0.02);
	a.setCoefficient(50.);
	if (a.getCoefficient() == 50.)
		std::cout << "\033[1;32msetCoefficient() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setCoefficient() does not perform as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto en = heatDissipationEngine<dsType, valType>(60.f, true, 0.001);
	auto l1 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 500, 1000, 1024));
	auto l2 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 1500, 1000, 1024));
	en.addAirLayer(l1);
	en.addAirLayer(l2);
	
	for (int i = 0; i < 10; ++i)
		en.step();

	int threads = omp_get_max_threads();
	std::cout << "\033[1;37mExecution time using " << threads << " threads for 10 steps of two grids with 1024 cell wide faces is: \033[1;33m" << en.getRunTimePassed() << " \033[1;37mseconds." << std::endl; 
	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'heatDissipationEngine'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
	testGetters();
	
	std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
	testSetters();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'heatDissipationEngine'\033[0m complete." << std::endl;
}