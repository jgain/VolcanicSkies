/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/
#include "airLayer.h"
#include <fstream>
#include <iostream>
#include "planet.h"
#include "terrain.h"
#include "moistureUptakeEngine.h"
#include <vector>

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
typedef PWM::PWMDataStructure::square2DArray<valType2> dsSType;
using namespace PWM::Engine;

void testConstructors(){
	auto a = moistureUptakeEngine<dsType, dsSType, valType, valType2>(60.f, true);
	if (a.getTerrain() == nullptr && a.getAirLayers().size() == 0)
		std::cout << "\033[1;32mDefault constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Default constructor does not perform as expected!\033[0m" << std::endl;


	auto p = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto t1 = PWM::Model::terrain<dsType, dsSType, valType, valType2>(p, 4);
	auto t = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t1);
	auto l = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(p, 100, 100, 1024));

	auto b = moistureUptakeEngine<dsType, dsSType, valType, valType2>(60.f, true);
	b.setTerrain(t);
	b.addAirLayer(l);
	if (b.getTerrain() == t && b.getAirLayers()[0] == l)
		std::cout << "\033[1;32mInitialising constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Initialising constructor does not perform as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	
	std::shared_ptr<PWM::Model::terrain<dsType, dsSType, valType, valType2>> t = nullptr;
	auto t1 = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, 1024);
	auto v = std::vector<std::pair<valType, bool>>();
	for (int i = 0; i < 1024 * 1024; ++i)
		v.push_back(std::make_pair(i / 10.f, i % 2));
	t1.init(v);
	t = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t1);

	auto a = moistureUptakeEngine<dsType, dsSType, valType, valType2>(60.f, true);
	a.setTerrain(t);

	auto l1 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 100, 100, 1024));
	auto l2 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 550, 1000, 1024));
	auto l3 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 1550, 1000, 1024));
	auto l4 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 2550, 1000, 1024));
	auto l5 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 3550, 1000, 1024));
	auto l6 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 4550, 1000, 1024));
	
	a.addAirLayer(l1);
	a.addAirLayer(l2);
	a.addAirLayer(l3);
	a.addAirLayer(l4);
	a.addAirLayer(l5);
	a.addAirLayer(l6);
	for (int i = 0; i < 10; ++i){
		a.step();
	}
	std::cout << "Time taken for 10 steps of moisture uptake for " << 1024 << "x" << 1024 << " grid: " << a.getRunTimePassed() << " seconds." << std::endl;
	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'moistureUptakeEngine'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'moistureUptakeEngine'\033[0m complete." << std::endl;
}