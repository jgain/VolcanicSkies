#include <iostream>
#include "position.h"
#include <tuple>

using namespace PWM::PWMDataStructure;
void testConstructors(){
    auto a = Position();
    bool p1, p2, p3;
    p1 = a.getX() == 0;
    p2 = a.getY() == 0;
    p3 = a.getZ() == 0;
    if (p1 && p2 && p3)
 	    std::cout << "\033[1;32mDefault constructor works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Default constructor does not perform as expected!\033[0m" << std::endl;

    auto b = Position(3, 4, 5);
    p1 = b.getX() == 3;
    p2 = b.getY() == 4;
    p3 = b.getZ() == 5;
    if (p1 && p2 && p3)
 	    std::cout << "\033[1;32mParametered constructor works as expected.\033[0m" << std::endl;
    else
	    std::cout << "\033[1;31mError! Parametered constructor does not perform as expected!\033[0m" << std::endl;
}

void testComparators(){
	auto a = Position(3, 4, 5);
    auto b = Position(3, 4, 5);
    if (a == b)
        std::cout << "\033[1;32mEquality operator works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator does not work as expected!\033[0m" << std::endl;
}

void testGetters(){
	auto a = Position(3.141, 2.786, 6.282);
    if(a.getX() == 3.141)
  	    std::cout << "\033[1;32mgetX() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getX() does not perform as expected!\033[0m" << std::endl;

    if(a.getY() == 2.786)
  	    std::cout << "\033[1;32mgetY() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getY() does not perform as expected!\033[0m" << std::endl;

    if(a.getZ() == 6.282)
  	    std::cout << "\033[1;32mgetZ() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getZ() does not perform as expected!\033[0m" << std::endl;

}

void testSetters(){
	auto a = Position(3, 4, 5);
    auto b = a;
    auto c = a;
    auto d = a;
    auto e = Position(8, 6, 4);

    a.move(e);
    if (a == e)
        std::cout << "\033[1;32mmove(Position& other) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! move(Position& other) does not perform as expected!\033[0m" << std::endl;

    double v1, v2, v3;
    v1 = 5;
    v2 = 2;
    v3 = -1;
    auto t1 = std::make_tuple(v1, v2, v3);
    b.move(t1);
    if (b == e)
        std::cout << "\033[1;32mmove(std::tuple<double, double, double>& distance) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! move(std::tuple<double, double, double>& distance) does not perform as expected!\033[0m" << std::endl;


    v1 = 10;
    v2 = 4;
    v3 = -2;
    t1 = std::make_tuple(v1, v2, v3);
    c.move(t1, 0.5);
    if (c == e)
        std::cout << "\033[1;32mmove(std::tuple<double, double, double>& velocity, double dt) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! move(std::tuple<double, double, double>& velocity, double dt) does not perform as expected!\033[0m" << std::endl;

    d.move(8, 6, 4);
    if (d == e)
        std::cout << "\033[1;32mmove(double newX, double newY, double newZ) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! move(double newX, double newY, double newZ) does not perform as expected!\033[0m" << std::endl;

}

void testMiscellany(){
	auto p0 = Position();
    auto a = Position(2, 3, 6);
    if (a.getDistanceTo(p0) == 7)
        std::cout << "\033[1;32mgetDistanceTo(Position& other) works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getDistanceTo(Position& other) does not perform as expected!\033[0m" << std::endl;

}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'Position'\033[0m:" << std::endl;
	
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
	
	std::cout << "\nTesting of class \033[1;34m'Position'\033[0m complete." << std::endl;
}