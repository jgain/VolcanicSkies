/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/
#include <iostream>
#include "planet.h"
#include "square2DArray.h"
#include <string>
#include "terrain.h"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
typedef PWM::PWMDataStructure::square2DArray<valType2> dsSType;
using namespace PWM::Model;

void testConstructors(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = terrain<dsType, dsSType, valType, valType2>(p, 128);
	bool p1, p2, p3;
	p1 = a.getPlanet() == p;
	p2 = a.getElevation().getWidth() == 128;
	p3 = a.getTemperature(1) == p->getAveTemp();
	if (p1 && p2 && p3)
		std::cout << "\033[1;32mConstructor(planet& p) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Constructor(planet& p) does not perform as expected!\033[0m" << std::endl;

	auto b = terrain<dsType, dsSType, valType, valType2>(p, 4);
	p1 = b.getPlanet() == p;
	p2 = b.getElevation().getWidth() == 4;
	p3 = b.getTemperature(1) == p->getAveTemp();
	if (p1 && p2 && p3)
		std::cout << "\033[1;32mConstructor(planet& p, int theta, int phi) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Constructor(planet& p, int theta, int phi) does not perform as expected!\033[0m" << std::endl;

	auto c(a);
	p1 = c.getPlanet() == p;
	p2 = c.getElevation().getWidth() == 128;
	p3 = c.getTemperature(1) == p->getAveTemp();
	if (p1 && p2 && p3)
		std::cout << "\033[1;32mCopy constructor(terrain& other) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Copy constructor(terrain& other) does not perform as expected!\033[0m" << std::endl;

	c = b;
	p1 = c.getPlanet() == p;
	p2 = c.getElevation().getWidth() == 4;
	p3 = c.getTemperature(1) == p->getAveTemp();
	if (p1 && p2 && p3)
		std::cout << "\033[1;32mCopy assignment(terrain& other) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Copy assignment(terrain& other) does not perform as expected!\033[0m" << std::endl;

}

void testComparators(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = terrain<dsType, dsSType, valType, valType2>(p, 4);
	auto b = terrain<dsType, dsSType, valType, valType2>(p, 4);
	auto c(a);
	if (a == a)
		std::cout << "\033[1;32mComparator correctly catches self comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch self comparison!\033[0m" << std::endl;
	
	if (a == b)
		std::cout << "\033[1;32mComparator correctly catches equal construction.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch equal construction!\033[0m" << std::endl;
	
	if (a == c)
		std::cout << "\033[1;32mComparator correctly catches copy construction.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch copy construction!\033[0m" << std::endl;

	valType t = 178;
	b.setTemperature(1, t);
	if (a != b)
		std::cout << "\033[1;32mComparator correctly catches different temperatures.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch different temperatures!\033[0m" << std::endl;
	
	b = a;
	t = 15;
	b.setMoisture(1, t);
	if (a != b)
		std::cout << "\033[1;32mComparator correctly catches different moistures.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch different moistures!\033[0m" << std::endl;
	
}

void testGetters(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = terrain<dsType, dsSType, valType, valType2>(p, 4);
	valType t = 178;
	
	if (a.getPlanet() == p)
		std::cout << "\033[1;32mgetPlanet() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getPlanet() does not work as expected!\033[0m" << std::endl;
	
	
	a.setTemperature(1, t);
	if (a.getTemperature().getData(1) == 178)
		std::cout << "\033[1;32mgetTemperature() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getTemperature() does not work as expected!\033[0m" << std::endl;
	
	if (a.getTemperature(1) == 178)
		std::cout << "\033[1;32mgetTemperature(int index) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getTemperature(int index) does not work as expected!\033[0m" << std::endl;
	
	t = 4;
	a.setMoisture(1, t);
	if (a.getMoisture().getData(1) == 4)
		std::cout << "\033[1;32mgetMoisture() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getMoisture() does not work as expected!\033[0m" << std::endl;
	
	if (a.getMoisture(1) == 4)
		std::cout << "\033[1;32mgetMoisture(int index) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getMoisture(int index) does not work as expected!\033[0m" << std::endl;
	

	valType2 t1 = "ICE";
	valType2 t2 = "VEGETATION";
	a.setTerrainType(1, t1);
	a.setTerrainType(5, t2);
	if (a.getTerrainType().getData(1) == "ICE")
		std::cout << "\033[1;32mgetTerrainType() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getTerrainType() does not work as expected!\033[0m" << std::endl;
	
	if (a.getTerrainType(5) == "VEGETATION")
		std::cout << "\033[1;32mgetTerrainType(int index) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getTerrainType(int index) does not work as expected!\033[0m" << std::endl;
	
}

void testSetters(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = terrain<dsType, dsSType, valType, valType2>(p, 4);
	valType t1 = 178;
	valType2 t2 = "ICE";
	a.setTemperature(1, t1);
	if (a.getTemperature(1) == 178)
		std::cout << "\033[1;32msetTemperature(int index, V& newTemp) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! setTemperature(int index, V& newTemp) does not work as expected!\033[0m" << std::endl;
	
	t1 = 89;
	a.setMoisture(1, t1);
	if (a.getMoisture(1) == 89)
		std::cout << "\033[1;32msetMoisture(int index, V& newTemp) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! setMoisture(int index, V& newTemp) does not work as expected!\033[0m" << std::endl;
	
	a.setTerrainType(1, t2);
	if (a.getTerrainType(1) == "ICE")
		std::cout << "\033[1;32msetTerrainType(int index, VV& newType) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! setTerrainType(int index, VV& newType) does not work as expected!\033[0m" << std::endl;

}

void testMiscellany(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = terrain<dsType, dsSType, valType, valType2>(p, 8);
	auto v = std::vector<std::pair<valType, bool>>();
	for (int i = 0; i < 64; ++i)
		v.push_back(std::make_pair(i, i % 2));
	a.init(v);

	if (a.getTerrainType(1) == "WATER" && a.getTerrainType(2) == "SOIL")
		std::cout << "\033[1;32minit() uses the bool value correctly.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! init() does not give expected terrain type values!\033[0m" << std::endl;
	
	if (a.getElevation(4) == 4)
		std::cout << "\033[1;32minit() gives expected elevation.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! init() does not give expected elevations!\033[0m" << std::endl;


	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'Terrain'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting comparators:\033[0m" << std::endl;
	testComparators();
	
	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
	testGetters();
	
	std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
	testSetters();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'Terrain'\033[0m complete." << std::endl;
}