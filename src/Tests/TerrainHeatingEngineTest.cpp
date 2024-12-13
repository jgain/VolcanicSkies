/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/
#include <fstream>
#include <iostream>
#include "planet.h"
#include "sun.h"
#include "sunEngine.h"
#include "terrain.h"
#include "terrainHeatingEngine.h"
#include <vector>

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
typedef PWM::PWMDataStructure::square2DArray<valType2> dsSType;
using namespace PWM::Engine;

void testConstructors(){
	auto a = terrainHeatingEngine<dsType, dsSType, valType, valType2>(60.f, true);
	if (a.getTerrain() == nullptr && a.getSunList().size() == 0)
		std::cout << "\033[1;32mDefault constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Default constructor does not perform as expected!\033[0m" << std::endl;


	auto p = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	auto t1 = PWM::Model::terrain<dsType, dsSType, valType, valType2>(p, 4);
	auto t = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t1);
	auto s = std::make_shared<PWM::Model::sun<valType>>(PWM::Model::sun<valType>("../resources/SunSol.json"));
	
	auto b = terrainHeatingEngine<dsType, dsSType, valType, valType2>(60.f, true, t);
	b.addSun(s);
	if (b.getTerrain() == t && b.getSunList()[0] == s)
		std::cout << "\033[1;32mInitialising constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Initialising constructor does not perform as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
	
	std::shared_ptr<PWM::Model::terrain<dsType, dsSType, valType, valType2>> t = nullptr;
	auto t1 = PWM::Model::terrain<dsType, dsSType, valType, valType2>(planetos, 256);
	auto v = std::vector<std::pair<valType, bool>>();
	for (int i = 0; i < 256 * 256; ++i)
		v.push_back(std::make_pair(i, i % 2));
	t1.init(v);
	t = std::make_shared<PWM::Model::terrain<dsType, dsSType, valType, valType2>>(t1);

	auto s = std::make_shared<PWM::Model::sun<valType>>(PWM::Model::sun<valType>("../resources/SunSol.json"));
	
	auto a = terrainHeatingEngine<dsType, dsSType, valType, valType2>(60.f, true, t);
	auto b = sunEngine<valType>(60.f, true);
	a.addSun(s);
	b.addSun(s);
	for (int i = 0; i < 10; ++i){
		a.step();
		b.step();
	}
	std::cout << "Time taken for 10 steps of groundHeating for " << 256 << "x" << 256 << " grid: " << a.getRunTimePassed() << " seconds." << std::endl;
	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'terrainHeatingEngine'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'terrainHeatingEngine'\033[0m complete." << std::endl;
}