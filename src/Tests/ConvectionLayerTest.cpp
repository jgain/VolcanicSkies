/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/

#include "convectionLayer.h"
#include <iostream>
#include "square2DArray.h"

typedef double valType;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
using namespace PWM::Model;

void testConstructors(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = convectionLayer<dsType, valType>(p, 200, 100, 8, 20000.f);
	bool p1, p2, p3, p4, p5;
	p1 = a.getPlanet() == p;
	p2 = a.getBottomHeight() == 100;
	p3 = a.getTopHeight() == 200;
	p4 = a.getVerticalVelocities().size() == 64;
	p5 = a.getLayerHeight() == 150;
	if (p1 && p2 && p3 && p4 && p5)
		std::cout << "\033[1;32mConstructor works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Constructor does not perform as expected!\033[0m" << std::endl;
}

void testComparators(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = convectionLayer<dsType, valType>(p, 200, 100, 4, 20000.f);
	auto b = convectionLayer<dsType, valType>(p, 200, 100, 4, 20000.f);
	
	if (a == a)
		std::cout << "\033[1;32mComparator correctly catches self comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch self comparison!\033[0m" << std::endl;
	
	if (a == b)
		std::cout << "\033[1;32mComparator correctly catches equal construction.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch equal construction!\033[0m" << std::endl;
	
	double t = 89;
	b.setVerticalVelocity(1, t);
	if (a != b)
		std::cout << "\033[1;32mComparator correctly catches different vertical velocities.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch different vertical velocities!\033[0m" << std::endl;
	
	auto c = convectionLayer<dsType, valType>(p, 200, 0, 4, 20000.f);
	if (!(a < a))
		std::cout << "\033[1;32mLess than operator correctly fails self comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Less than operator passes self comparison!\033[0m" << std::endl;
	
	if (!(a < b))
		std::cout << "\033[1;32mLess than operator correctly fails equal comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Less than operator passes equal comparison!\033[0m" << std::endl;
	
	if (c < a)
		std::cout << "\033[1;32mLess than operator correctly passes comparison of lhs < rhs.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Less than operator fails comparison of lhs < rhs!\033[0m" << std::endl;
	
}

void testGetters(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = convectionLayer<dsType, valType>(p, 200, 100, 4, 20000.f);
	if (a.getPlanet() == p)
		std::cout << "\033[1;32mgetPlanet() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getPlanet() does not work as expected!\033[0m" << std::endl;
	
	if (a.getBottomHeight() == 100)
		std::cout << "\033[1;32mgetBottomHeight() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getBottomHeight() does not work as expected!\033[0m" << std::endl;
	
	if (a.getTopHeight() == 200)
		std::cout << "\033[1;32mgetTopHeight() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getTopHeight() does not work as expected!\033[0m" << std::endl;
	
	double t = 89;
	a.setVerticalVelocity(9, t);
	if (a.getVerticalVelocities().getData(9) == 89)
		std::cout << "\033[1;32mgetVerticalVelocities() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getVerticalVelocities() does not work as expected!\033[0m" << std::endl;
	
	if (a.getVerticalVelocity(9) == 89)
		std::cout << "\033[1;32mgetVerticalVelocity(const size_t index) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getVerticalVelocity(const size_t index) does not work as expected!\033[0m" << std::endl;
	
	if (a.getVerticalVelocity(2, 1) == 89)
		std::cout << "\033[1;32mgetVerticalVelocity(const size_t face, const size_t i, const size_t j) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getVerticalVelocity(const size_t face, const size_t i, const size_t j) does not work as expected!\033[0m" << std::endl;
	
	if (a.getLayerHeight() == 150)
		if (!(a < a))
		std::cout << "\033[1;32mgetLayerHeight() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getLayerHeight() does not work as expected!\033[0m" << std::endl;
	
}

void testSetters(){
	auto p = std::make_shared<planet>(planet("../resources/PlanetEarth.json"));
	auto a = convectionLayer<dsType, valType>(p, 200, 100, 4, 20000.f);
	double t = 89;
	a.setVerticalVelocity(1, t);
	if (a.getVerticalVelocity(1) == 89)
		std::cout << "\033[1;32msetVerticalVelocity(int index, V newVal) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! setVerticalVelocity(int index, V newVal) does not work as expected!\033[0m" << std::endl;
	
	a.setVerticalVelocity(3, 1, t);
	if (a.getVerticalVelocity(3, 1) == 89)
		std::cout << "\033[1;32msetVerticalVelocity(const size_t face, const size_t i, const size_t j, V newVal) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! setVerticalVelocity(const size_t face, const size_t i, const size_t j, V newVal) does not work as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'ConvectionLayer'\033[0m:" << std::endl;
	
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
	
	std::cout << "\nTesting of class \033[1;34m'ConvectionLayer'\033[0m complete." << std::endl;
}