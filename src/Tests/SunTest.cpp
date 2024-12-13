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
#include <iostream>
#include "sun.h"

using namespace PWM::Model;

void testConstructors(){
	auto a = sun<double>();
	if (a.getApparentDeclination() == 0 && a.getApparentRightAscension() == 0 && a.getSeasonDirection() == 1)
		std::cout << "\033[1;32mBlank constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Blank constructor does not perform as expected!\033[0m" << std::endl;

	auto b = sun<double>("../resources/SunSol.json");
	bool p1, p2, p3, p4, p5, p6, p7;
	p1 = b.getPower() == 1350;
	p2 = b.getColour() == 550;
	p3 = b.getDistance() == 1;
	p4 = b.getApparentDeclination() == 0 && b.getApparentRightAscension() == 0 && a.getSeasonDirection() == 1;
	p5 = b.getAxialTilt() == 23.5;
	p6 = b.getTropicalYearLength() == 31556923.488;
	p7 = b.getSolarDayLength() == 86400;
	if (p1 && p2 && p3 && p4 && p5 && p6 && p7)
		std::cout << "\033[1;32mJSON constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! JSON constructor does not perform as expected!\033[0m" << std::endl;
}

void testComparators(){
	auto a = sun<double>("../resources/SunSol.json");
	auto b = sun<double>("../resources/SunSol.json");

	if (a == a)
		std::cout << "\033[1;32mEquality operator works works as expected with self comparison.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Self comparison does not work as expected!\033[0m" << std::endl;

	if (a == b)
		std::cout << "\033[1;32mEquality operator works as expected for equal construction.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Equal construction does not give positive result from equality operator!\033[0m" << std::endl;
}

void testGetters(){
	auto a = sun<double>("../resources/SunSol.json");
	
	if (a.getPower() == 1350)
		std::cout << "\033[1;32mgetPower() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getPower() does not perform as expected!\033[0m" << std::endl;

	if (a.getColour() == 550)
		std::cout << "\033[1;32mgetColour() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getColour() does not perform as expected!\033[0m" << std::endl;

	if (a.getDistance() == 1)
		std::cout << "\033[1;32mgetDistance() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getDistance() does not perform as expected!\033[0m" << std::endl;
	
	if (a.getApparentDeclination() == 0)
		std::cout << "\033[1;32mgetApparentDeclination() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getApparentDeclination() does not perform as expected!\033[0m" << std::endl;
	
	if (a.getSeasonDirection() == 1)
		std::cout << "\033[1;32mgetSeasonDirection() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getSeasonDirection() does not perform as expected!\033[0m" << std::endl;

	if (a.getApparentRightAscension() == 0)
		std::cout << "\033[1;32mgetApparentRightAscension() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getApparentRightAscension() does not perform as expected!\033[0m" << std::endl;
	
	if (a.getAxialTilt() == 23.5)
		std::cout << "\033[1;32mgetAxialTilt() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getAxialTilt() does not perform as expected!\033[0m" << std::endl;
	
	if (a.getTropicalYearLength() == 31556923.488)
		std::cout << "\033[1;32mgetTropicalYearLength() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getTropicalYearLength() does not perform as expected!\033[0m" << std::endl;
	
	if (a.getSolarDayLength() == 86400)
		std::cout << "\033[1;32mgetSolarDayLength() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getSolarDayLength() does not perform as expected!\033[0m" << std::endl;
}

void testSetters(){
	auto a = sun<double>("../resources/SunSol.json");
	
	a.setPower(1500);
	if (a.getPower() == 1500)
		std::cout << "\033[1;32msetPower(T pow) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setPower(T pow) does not perform as expected!\033[0m" << std::endl;

	a.setColour(750);
	if (a.getColour() == 750)
		std::cout << "\033[1;32msetColour(T col) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setColour(T col) does not perform as expected!\033[0m" << std::endl;

	a.setDistance(2);
	if (a.getDistance() == 2)
		std::cout << "\033[1;32msetDistance(T dis) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setDistance(T dis) does not perform as expected!\033[0m" << std::endl;
	
	a.setApparentDeclination(0.002);
	if (a.getApparentDeclination() == 0.002)
		std::cout << "\033[1;32msetApparentDeclination(T dec) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setApparentDeclination(T dec) does not perform as expected!\033[0m" << std::endl;
	
	a.flipSeasonDirection();
	if (a.getSeasonDirection() == -1)
		std::cout << "\033[1;32msetSeasonDirection() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setSeasonDirection() does not perform as expected!\033[0m" << std::endl;

	a.setApparentRightAscension(0.01);
	if (a.getApparentRightAscension() == 0.01)
		std::cout << "\033[1;32msetApparentRightAscension(T ra) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setApparentRightAscension(T ra) does not perform as expected!\033[0m" << std::endl;
	
	a.setAxialTilt(25.19);
	if (a.getAxialTilt() == 25.19)
		std::cout << "\033[1;32msetAxialTilt(T axt) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setAxialTilt(T axt) does not perform as expected!\033[0m" << std::endl;
	
	a.setTropicalYearLength(59355072.);
	if (a.getTropicalYearLength() == 59355072.)
		std::cout << "\033[1;32msetTropicalYearLength() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setTropicalYearLength() does not perform as expected!\033[0m" << std::endl;
	
	a.setSolarDayLength(88776);
	if (a.getSolarDayLength() == 88776)
		std::cout << "\033[1;32msetSolarDayLength() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setSolarDayLength() does not perform as expected!\033[0m" << std::endl;
}

void testMiscellany(){
	
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'sun'\033[0m:" << std::endl;
	
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
	
	std::cout << "\nTesting of class \033[1;34m'sun'\033[0m complete." << std::endl;
}