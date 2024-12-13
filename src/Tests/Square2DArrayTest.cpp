/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
 */
#include <iostream>
#include "square2DArray.h"

using namespace PWM::PWMDataStructure;

void testConstructors(){
    auto a = square2DArray<double>(4);
    if (a.getWidth() == 4)
       std::cout << "\033[1;32mBasic constructor appears to work as expected.\033[0m" << std::endl;
    else
       std::cout << "\033[1;31mError! Basic constructor does not behave as expected!\033[0m" << std::endl;

    if (a.getData(1, 1) == 0 && a.getSize() == 50000.f)
        std::cout << "\033[1;32mInitialisation works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Initialisation does not perform as expected!\033[0m" << std::endl;

    auto b = square2DArray<double>(8, 20000.f);
    if (b.getWidth() == 8 && b.getData(10) == 0 && b.getSize() == 20000.f)
        std::cout << "\033[1;32mSized constructor appears to work as expected.\033[0m" << std::endl;
    else
       std::cout << "\033[1;31mError! Sized constructor does not behave as expected!\033[0m" << std::endl;
}

void testComparators(){
	auto a = square2DArray<double>(4);
    auto b = square2DArray<double>(4);
    auto c = square2DArray<double>(4, 20000.f);
    if (a == a)
		std::cout << "\033[1;32mComparator correctly catches self comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch self comparison!\033[0m" << std::endl;
	
	if (a == b)
		std::cout << "\033[1;32mComparator correctly catches equal construction.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch equal construction!\033[0m" << std::endl;
	
    if (a != c)
        std::cout << "\033[1;32mComparator correctly catches different sizes.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch different sizes!\033[0m" << std::endl;
	
	a.setData(1, 1, 5);
    if (a != b)
		std::cout << "\033[1;32mComparator correctly catches different data values.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Comparator fails to catch different data values!\033[0m" << std::endl;
}

void testGetters(){
	auto a = square2DArray<double>(4, 20000);
    a.setData(1, 1, 5);
    if (a.getData(1, 1) == 5)
        std::cout << "\033[1;32mget(const size_t i, const size_t j) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! get(const size_t i, const size_t j) does not perform as expected!\033[0m" << std::endl;
    
    a.setData(1, 1, 1);
    a.setData(1, 2, 2);
    a.setData(2, 1, 3);
    a.setData(2, 2, 4);
    if (a.getInterpolated(1.5f, 1.5f) == 2.5f)
        std::cout << "\033[1;32mgetInterpolated(const float x, const float y) works as expected.\033[0m" << std::endl;
    else{
        std::cout << "\033[1;31mError! getInterpolated(const float x, const float y) does not perform as expected!\033[0m" << std::endl;
        std::cout << "Expected: " << 2.5f << std::endl;
        std::cout << "Actual: " << a.getInterpolated(1.5f, 1.5f) << std::endl;
    }

    if (a.getWidth() == 4)
        std::cout << "\033[1;32mgetWidth() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getWidth() does not work as expected!\033[0m" << std::endl;

    if (a.getSize() == 20000.f)
        std::cout << "\033[1;32mgetSize() works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! getSize() does not work as expected!\033[0m" << std::endl;

    double x = 20000.f / 4;
    if (a.gridLength() == x)
        std::cout << "\033[1;32mgridLength() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! gridLength() does not perform as expected!\033[0m" << std::endl;

    auto t = Coordinate(0, 0);
    if (a.getCoordinates() == t)
        std::cout << "\033[1;32mgetCoordinates() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getCoordinates() does not perform as expected!\033[0m" << std::endl;
}

void testSetters(){
	auto a = square2DArray<double>(4);
    a.setData(1, 1, 5);
    if (a.getData(1, 1) == 5)
        std::cout << "\033[1;32mset(const size_t i, const size_t j, const T& val) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! set(const size_t i, const size_t j, const T& val) does not perform as expected!\033[0m" << std::endl;
    
    auto t = Coordinate(46.196, -122.184);
    a.setCoord(t);
    auto x = Coordinate(46.196, -122.184);
    if (a.getCoordinates() == x)
        std::cout << "\033[1;32msetCoord(const Coordinate& nCoord) works as expected.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! setCoord(const Coordinate& nCoord) does not work as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	auto a = square2DArray<double>(4, 50000.f);
    auto b = square2DArray<double>(4, 50000.f);
    a.randomInit(85000, 108000);
    b.randomInit(85000, 108000);
    auto c = square2DArray<double>(4, 50000.f);
    a.laplacian(c);
    std::cout << "\033[1;32mLaplacian(square2DArray<T>& res) performed apparently without issue.\033[0m" << std::endl;

    a.gaussSeidelSmoothBS(c, b, 50000.f);
    std::cout << "\033[1;32mgaussSeidelSmoothBS(square2DArray<T>& res, const square2DArray<T>& other, const T gridSize) performed apparently without issue.\033[0m" << std::endl;

    b.gaussSeidelSmoothBS(c, 10, 50000.f);
    std::cout << "\033[1;32mgaussSeidelSmoothBS(const square2DArray<T>& other, const int nb, const T gridSize) performed apparently without issue.\033[0m" << std::endl;
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'square2DArray'\033[0m:" << std::endl;
	
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
	
	std::cout << "\nTesting of class \033[1;34m'square2DArray'\033[0m complete." << std::endl;
}