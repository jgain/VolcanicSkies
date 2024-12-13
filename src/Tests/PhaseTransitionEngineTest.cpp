/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/

#include "airLayer.h"
#include <iostream>
#include <omp.h>
#include "phaseTransitionEngine.h"
#include "planet.h"
#include "square2DArray.h"

typedef double valType;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
using namespace PWM::Engine;

void testConstructors(){
	auto a = phaseTransitionEngine<dsType, valType>(60.f, true);
	if (a.getAirLayers().size() == 0)
		std::cout << "\033[1;32mBasic constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Basic constructor does not perform as expected!\033[0m" << std::endl;

	auto b = phaseTransitionEngine<dsType, valType>(60.f, true, true, true, true, true, 0.001, 0.001);
	bool p1 = b.getDoEvaporation() && b.getDoCondensation();
	bool p2 = b.getDoEvaporationHX() && b.getDoCondensationHX();
	bool p3 = b.getCondensationCoefficient() == 0.001 && b.getEvaporationCoefficient() == 0.001;
	bool p4 = b.getAirLayers().size() == 0;
	if (p1 && p2 && p3 && p4)
		std::cout << "\033[1;32mFull constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Full constructor does not perform as expected, " << p1 << ", " << p2 << ", " << p3 << ", and " << p4 << "!\033[0m" << std::endl;
}

void testGetters(){
	auto a = phaseTransitionEngine<dsType, valType>(60.f, true, true, false, false, true, 0.4f, 0.1f);
	if (a.getDoEvaporation())
		std::cout << "\033[1;32mgetDoEvaporation() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getDoEvaporation() does not perform as expected!\033[0m" << std::endl;

	if (!a.getDoEvaporationHX())
		std::cout << "\033[1;32mgetDoEvaporationHX() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getDoEvaporationHX() does not perform as expected!\033[0m" << std::endl;

	if (!a.getDoCondensation())
		std::cout << "\033[1;32mgetDoCondensation() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getDoCondensation() does not perform as expected!\033[0m" << std::endl;

	if (a.getDoCondensationHX())
		std::cout << "\033[1;32mgetDoCondensationHX() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getDoCondensationHX() does not perform as expected!\033[0m" << std::endl;

	if (a.getEvaporationCoefficient() == 0.4f)
		std::cout << "\033[1;32mgetEvaporationCoefficient() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getEvaporationCoefficient() does not perform as expected!\033[0m" << std::endl;

	if (a.getCondensationCoefficient() == 0.1f)
		std::cout << "\033[1;32mgetCondensationCoefficient() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getCondensationCoefficient() does not perform as expected!\033[0m" << std::endl;
}

void testSetters(){
	auto a = phaseTransitionEngine<dsType, valType>(60.f, true, true, false, false, true, 0.4f, 0.1f);
	a.setDoEvaporation(false);
	if (!a.getDoEvaporation())
		std::cout << "\033[1;32msetDoEvaporation() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setDoEvaporation() does not perform as expected!\033[0m" << std::endl;

	a.setDoEvaporationHX(true);
	if (a.getDoEvaporationHX())
		std::cout << "\033[1;32msetDoEvaporationHX() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setDoEvaporationHX() does not perform as expected!\033[0m" << std::endl;

	a.setDoCondensation(true);
	if (a.getDoCondensation())
		std::cout << "\033[1;32msetDoCondensation() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setDoCondensation() does not perform as expected!\033[0m" << std::endl;

	a.setDoCondensationHX(false);
	if (!a.getDoCondensationHX())
		std::cout << "\033[1;32msetDoCondensationHX() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setDoCondensationHX() does not perform as expected!\033[0m" << std::endl;

	a.setEvaporationCoefficient(24.);
	if (a.getEvaporationCoefficient() == 24.)
		std::cout << "\033[1;32msetEvaporationCoefficient() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setEvaporationCoefficient() does not perform as expected!\033[0m" << std::endl;

	a.setCondensationCoefficient(69.);
	if (a.getCondensationCoefficient() == 69.)
		std::cout << "\033[1;32msetCondensationCoefficient() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setCondensationCoefficient() does not perform as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto en = phaseTransitionEngine<dsType, valType>(60.f, true, true, true, true, true, 0.001, 0.001);
	auto l1 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 500, 1000, 1024));
	auto l2 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 1500, 1000, 1024));
	l1->randomInit();
	l2->randomInit();
	en.addAirLayer(l1);
	en.addAirLayer(l2);
	
	for (int i = 0; i < 10; ++i)
		en.step();

	int threads = omp_get_max_threads();
	std::cout << "\033[1;37mExecution time using " << threads << " threads for 10 steps of two grids with 1024 cell wide faces is: \033[1;33m" << en.getRunTimePassed() << " \033[1;37mseconds." << std::endl; 
	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'PhaseTransitionEngine'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
	testGetters();
	
	std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
	testSetters();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'PhaseTransitionEngine'\033[0m complete." << std::endl;
}