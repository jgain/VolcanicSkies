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
// File to run tests for the Coordinate class.
// Cilliers Pretorius
// 14 May 2021
#include <iostream>
#include <utility>
#include "Coordinate.h"

using namespace PWM::PWMDataStructure;
void testConstructors(){
	Coordinate a = Coordinate(-31, 18);
	bool p1, p2;
	p1 = a.getLatitude() == -31;
	p2 = a.getLongitude() == 18;
	if (p1 && p2)
		std::cout << "\033[1;32mBasic constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Basic constructor does not perform as expected!\033[0m" << std::endl;

	bool p3;
	Coordinate b(a);
	p1 = a.getLatitude() == b.getLatitude();
	p2 = a.getLongitude() == b.getLongitude();
	p3 = a == b;
	if (p1 && p2 && p3)
		std::cout << "\033[1;32mCopy construction works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Copy construction is broken!\033[0m" << std::endl;

	Coordinate c = a;
	p1 = a.getLatitude() == c.getLatitude();
	p2 = a.getLongitude() == c.getLongitude();
	p3 = a == c;
	if (p1 && p2 && p3)
		std::cout << "\033[1;32mCopy assignment works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Copy assignment is broken!\033[0m" << std::endl;
}

void testSetters(){
	Coordinate a = Coordinate(-31, 18);
	bool p1, p2;
	a.setLatitude(25.5);
	a.setLongitude(28.05);
	p1 = a.getLatitude() == 25.5;
	p2 = a.getLongitude() == 28.05;
	if (p1 && p2)
		std::cout << "\033[1;32mSetters work as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Setters do not perform as expected!\033[0m" << std::endl;
}

void testGetters(){
	Coordinate a = Coordinate(-31, 18);
	bool p1, p2;
	p1 = a.getLatitude() == -31;
	p2 = a.getLongitude() == 18;
	if (p1 && p2)
		std::cout << "\033[1;32mGetters work as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Getters do not perform as expected!\033[0m" << std::endl;
}

void testComparators(){
	Coordinate a = Coordinate(-31, 18);
	Coordinate b(a);
	if (a == b)
		std::cout << "\033[1;32mThe '==' operator works as expected for copied objects.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! The '==' operator is broken!\033[0m" << std::endl;
	
	Coordinate c = Coordinate(-31, 18);
	if (a == c)
		std::cout << "\033[1;32mThe '==' operator works as expected for new objects.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! The '==' operator is broken!\033[0m" << std::endl;

	Coordinate d = Coordinate(-31.002, 18);
	if (a != d)
		std::cout << "\033[1;32mThe '!=' operator works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! The '!=' operator is broken!\033[0m" << std::endl;

	Coordinate e = Coordinate(-31.0001, 18.0001);
	bool p1 = a != e;
	bool p2 = a.approxEquals(e);
	bool p3 = !(a.approxEquals(d));
	if (p1 && p2)
		std::cout << "\033[1;32mapproxEquals(const Coordinate& other) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! approxEquals(const Coordinate& other) is broken!\033[0m" << std::endl;

	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'Coordinate'\033[0m:\n" << std::endl;

	std::cout << "\033[1;37mTesting comparators:\033[0m" << std::endl;
	testComparators();

	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();

	std::cout << "\n\033[1;37mTesting getters:\033[0m" << std::endl;
	testGetters();

	std::cout << "\n\033[1;37mTesting setters:\033[0m" << std::endl;
	testSetters();

	std::cout << "\nTesting of class \033[1;34m'Coordinate'\033[0m complete." << std::endl;
}