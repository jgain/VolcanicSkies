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
#include "airLayer.h"
#include "square2DArray.h"
#include <iostream>

typedef double valType;
typedef PWM::PWMDataStructure::square2DArray<valType> dsType;
using namespace PWM::Model;
void testConstructors(){
	auto planetos = std::make_shared<planet>(planet("../resources//PlanetEarth.json"));
	auto a = airLayer<dsType, valType>(planetos, 100, 100, 4, 20000.);
	bool p1, p2, p3, p4;
	p1 = (a.getHeight() == 100.);
	p2 = (a.getThickness() == 100.);
	p3 = (a.getPlanet() == planetos);
	p4 = ((a.getObstacles().size() == 16) && (a.getVelocityTheta().size() == 16) && (a.getVelocityPhi().size() == 16) && (a.getTemperature().size() == 16) && (a.getMoisture().size() == 16) && (a.getCondensedWater().size() == 16) && (a.getPressure().size() == 16));
	if (p1 && p2 && p3 && p4)
		std::cout << "\033[1;32mDefault constructor works as expected.\033[0m" << std::endl;
	else if (!p1)
		std::cout << "\033[1;31mError! Default constructor gives wrong height!\033[0m" << std::endl;
	else if (!p2)
		std::cout << "\033[1;31mError! Default constructor gives wrong thickness!\033[0m" << std::endl;
	else if (!p3)
		std::cout << "\033[1;31mError! Default constructor gives wrong planet!\033[0m" << std::endl;
	else if (!p4)
		std::cout << "\033[1;31mError! Default constructor creates wrong data structures!\033[0m" << std::endl;
	
	auto b(a);
	if (a == b)
		std::cout << "\033[1;32mCopy constructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Copy constructor does not perform as expected!\033[0m" << std::endl;
	b.setHeight(1);
	b = a;
	if (a == b)
		std::cout << "\033[1;32mCopy assignment works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Copy assignment does not perform as expected!\033[0m" << std::endl;
	
	auto d = std::make_shared<airLayer<dsType, valType>>(a);
	auto c(d);
	double pi = 3.1415926535;
	auto x = std::make_shared<dsType>(dsType(4));
	auto y = std::make_shared<dsType>(dsType(4));
	x->setData(0, pi);
	c->swapVels(x, y);
	if (a != *c)
		std::cout << "\033[1;32mCopy constructor works as expected and does not create a new pointer.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Copy constructor does not perform as expected and creates new pointers!\033[0m" << std::endl;
}

void testComparators(){
	auto planetos = std::make_shared<planet>(planet("../resources//PlanetEarth.json"));
	auto a = airLayer<dsType, valType>(planetos, 100, 100, 4, 20000.);
	auto b = airLayer<dsType, valType>(planetos, 100, 100, 4, 20000.);
	if (a == a)
		std::cout << "\033[1;32mOperator== accepts self comparison as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Operator== rejects self comparison!\033[0m" << std::endl;
	
	if (a == b)
		std::cout << "\033[1;32mOperator== accepts equal construction as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Operator== rejects equal construction!\033[0m" << std::endl;
	
	auto c = a;
	c.setHeight(1);
	if (a != c)
		std::cout << "\033[1;32mOperator!= accepts height difference as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Operator!= rejects height difference!\033[0m" << std::endl;
	
	c = a;
	c.setThickness(1);
	if (a != c)
		std::cout << "\033[1;32mOperator!= accepts thickness difference as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Operator!= rejects thickness difference!\033[0m" << std::endl;

	c = a;
	if (!(a < a))
		std::cout << "\033[1;32mLess than operator correctly fails self comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Less than operator passes self comparison!\033[0m" << std::endl;
	
	if (!(a < b))
		std::cout << "\033[1;32mLess than operator correctly fails equal comparison.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Less than operator passes equal comparison!\033[0m" << std::endl;
	
	c.setHeight(1);
	if (c < a)
		std::cout << "\033[1;32mLess than operator correctly passes comparison of lhs < rhs.\033[0m" << std::endl;
    else
		std::cout << "\033[1;31mError! Less than operator fails comparison of lhs < rhs!\033[0m" << std::endl;
	
}

void testGetters(){
	auto planetos = std::make_shared<planet>(planet("../resources//PlanetEarth.json"));
	auto a = airLayer<dsType, valType>(planetos, 100, 100, 4, 20000.);
	if (a.getRadius() == 6371100.)
		std::cout << "\033[1;32mgetRadius() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getRadius() gives wrong result!\033[0m" << std::endl;

	if(a.getHeight() == 100.)
		std::cout << "\033[1;32mgetHeight() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getHeight() gives wrong result!\033[0m" << std::endl;

	if(a.getThickness() == 100.)
		std::cout << "\033[1;32mgetThickness() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getThickness() gives wrong result!\033[0m" << std::endl;
	
	if(a.getPlanet() == planetos)
		std::cout << "\033[1;32mgetPlanet() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getPlanet() gives wrong result!\033[0m" << std::endl;

	a.setTemperature(10.);
	if(a.getTemperature().getData(1) == 10.)
		std::cout << "\033[1;32mgetTemperature() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getTemperature() gives wrong result!\033[0m" << std::endl;

	if(a.getTemperature(1) == 10.)
		std::cout << "\033[1;32mgetTemperature(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getTemperature(const size_t index) gives wrong result!\033[0m" << std::endl;
	
	if(a.getMeanTemperature() == 10.)
		std::cout << "\033[1;32mgetMeanTemperature() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getMeanTemperature() gives wrong result!\033[0m" << std::endl;

	a.setMoisture(19.);
	if(a.getMoisture().getData(10) == 19.)
		std::cout << "\033[1;32mgetMoisture() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getMoisture() gives wrong result!\033[0m" << std::endl;

	if(a.getMoisture(10) == 19.)
		std::cout << "\033[1;32mgetMoisture(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getMoisture(const size_t index) gives wrong result!\033[0m" << std::endl;
	
	if(a.getCumulativeMoisture() == (19 * a.getMoisture().size()))
		std::cout << "\033[1;32mgetCumulativeMoisture() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getCumulativeMoisture() gives wrong result!\033[0m" << std::endl;

	auto t = dsType(4);
	auto v = dsType(4);
	double t1 = 1.;
	double t2 = 3.45;
	t.setData(1, t1);
	t.setData(2, 1, t2);
	v.setData(1, t1);
	a.setObstacles(t);
	if(a.getObstacles().getData(1) == 1.)
		std::cout << "\033[1;32mgetObstacles() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getObstacles() gives wrong result!\033[0m" << std::endl;

	if(a.getObstacles(1) == 1.)
		std::cout << "\033[1;32mgetObstacles(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getObstacles(const size_t index) gives wrong result!\033[0m" << std::endl;
	
	a.setVelocityTheta(v);
	if(a.getVelocityTheta().getData(1) == 1.)
		std::cout << "\033[1;32mgetVelocityTheta() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getVelocityTheta() gives wrong result!\033[0m" << std::endl;

	if(a.getVelocityTheta(1) == 1.)
		std::cout << "\033[1;32mgetVelocityTheta(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getVelocityTheta(const size_t index) gives wrong result!\033[0m" << std::endl;
	
	a.setVelocityPhi(v);
	if(a.getVelocityPhi().getData(1) == 1.)
		std::cout << "\033[1;32mgetVelocityPhi() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getVelocityPhi() gives wrong result!\033[0m" << std::endl;

	if(a.getVelocityPhi(1) == 1.)
		std::cout << "\033[1;32mgetVelocityPhi(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getVelocityPhi(const size_t index) gives wrong result!\033[0m" << std::endl;
	
	a.setCondensedWater(t);
	if(a.getCondensedWater().getData(1) == 1.)
		std::cout << "\033[1;32mgetCondensedWater() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getCondensedWater() gives wrong result!\033[0m" << std::endl;

	if(a.getCondensedWater(1) == 1.)
		std::cout << "\033[1;32mgetCondensedWater(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getCondensedWater(const size_t index) gives wrong result!\033[0m" << std::endl;
	
	a.setPressure(t);
	if(a.getPressure().getData(1) == 1.)
		std::cout << "\033[1;32mgetPressure() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getPressure() gives wrong result!\033[0m" << std::endl;

	if(a.getPressure(1) == 1.)
		std::cout << "\033[1;32mgetPressure(const size_t index) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! getPressure(const size_t index) gives wrong result!\033[0m" << std::endl;
}

void testSetters(){
	auto planetos = std::make_shared<planet>(planet("../resources//PlanetEarth.json"));
	auto a = airLayer<dsType, valType>(planetos, 100, 100, 4, 20000.f);
	a.setHeight(2.);
	if(a.getHeight() == 2.)
		std::cout << "\033[1;32msetHeight() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setHeight() gives wrong result!\033[0m" << std::endl;

	a.setThickness(3.);
	if(a.getThickness() == 3.)
		std::cout << "\033[1;32msetThickness() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setThickness() gives wrong result!\033[0m" << std::endl;

	a.setTemperature(10.);
	if(a.getTemperature().getData(1) == 10.)
		std::cout << "\033[1;32msetTemperature(double heat) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setTemperature(double heat) gives wrong result!\033[0m" << std::endl;

	a.setMoisture(19.);
	if(a.getMoisture().getData(10) == 19.)
		std::cout << "\033[1;32msetMoisture(double vap) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setMoisture(double vap) gives wrong result!\033[0m" << std::endl;

	auto t = dsType(4);
	auto v = dsType(4);
	double t1 = 1.;
	double t2 = 8.34;
	t.setData(1, t1);
	v.setData(1, t1);
	a.setObstacles(t);
	if(a.getObstacles().getData(1) == 1.)
		std::cout << "\033[1;32msetObstacles() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setObstacles() gives wrong result!\033[0m" << std::endl;

	a.setObstacles(1, t2);
	if(a.getObstacles(1) == t2)
		std::cout << "\033[1;32msetObstacles(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setObstacles(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;

	a.setVelocityTheta(v);
	if(a.getVelocityTheta().getData(1) == 1.)
		std::cout << "\033[1;32msetVelocityTheta() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setVelocityTheta() gives wrong result!\033[0m" << std::endl;

	a.setVelocityTheta(1, t2);
	if(a.getVelocityTheta(1) == t2)
		std::cout << "\033[1;32msetVelocityTheta(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setVelocityTheta(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;

	v.setData(1, t1);
	a.setVelocityPhi(v);
	if(a.getVelocityPhi().getData(1) == 1.)
		std::cout << "\033[1;32msetVelocityPhi() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setVelocityPhi() gives wrong result!\033[0m" << std::endl;

	a.setVelocityPhi(1, t2);
	if(a.getVelocityPhi(1) == t2)
		std::cout << "\033[1;32msetVelocityPhi(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setVelocityPhi(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;

	t.setData(1, t1);
	a.setTemperature(t);
	if(a.getTemperature().getData(1) == 1.)
		std::cout << "\033[1;32msetTemperature(T& temps) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setTemperature(T& temps) gives wrong result!\033[0m" << std::endl;

	a.setTemperature(1, t2);
	if(a.getTemperature(1) == t2)
		std::cout << "\033[1;32msetTemperature(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setTemperature(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;

	t.setData(1, t1);
	a.setMoisture(t);
	if(a.getMoisture().getData(1) == 1.)
		std::cout << "\033[1;32msetMoisture(T& moist) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setMoisture(T& moist) gives wrong result!\033[0m" << std::endl;

	a.setMoisture(1, t2);
	if(a.getMoisture(1) == t2)
		std::cout << "\033[1;32msetMoisture(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setMoisture(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;

	t.setData(1, t1);
	a.setCondensedWater(t);
	if(a.getCondensedWater().getData(1) == 1.)
		std::cout << "\033[1;32msetCondensedWater() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setCondensedWater() gives wrong result!\033[0m" << std::endl;

	a.setCondensedWater(1, t2);
	if(a.getCondensedWater(1) == t2)
		std::cout << "\033[1;32msetCondensedWater(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setCondensedWateryTheta(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;

	t.setData(1, t1);
	a.setPressure(t);
	if(a.getPressure().getData(1) == 1.)
		std::cout << "\033[1;32msetPressure() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setPressure() gives wrong result!\033[0m" << std::endl;

	a.setPressure(1, t2);
	if(a.getPressure(1) == t2)
		std::cout << "\033[1;32msetPressure(const size_t index, const V& val) works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! setPressure(const size_t index, const V& val) gives wrong result!\033[0m" << std::endl;
}

void testMiscellany(){
	auto planetos = std::make_shared<planet>(planet("../resources//PlanetEarth.json"));
	auto a = airLayer<dsType, valType>(planetos, 100, 100,4, 20000.f);
	auto t = dsType(4);
	auto t1 = dsType(4);
	double f = 8.8;
	t1 = 1.1;
	a.setPressure(t1);
	auto w = std::make_shared<dsType>(dsType(4));
	auto x = std::make_shared<dsType>(dsType(4));
	auto y = std::make_shared<dsType>(dsType(4));
	auto z = std::make_shared<dsType>(dsType(4));
	auto vec = std::vector<std::shared_ptr<dsType>>();
	vec.push_back(w);
	vec.push_back(x);
	vec.push_back(y);
	vec.push_back(z);
	(*z) = f;
	a.swapScalars(vec);
	if(a.getPressure().getData(0) == f)
		std::cout << "\033[1;32mswapScalars() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! swapScalars() gives wrong result!\033[0m" << std::endl;

	//Testing the random initialisation
	auto b = airLayer<dsType, valType>(planetos, 100, 100, 4, 20000.f);
	bool inRange = true;
	std::pair<double, double> vMinMax = std::make_pair(-planetos->getMaxWindSpeed(), planetos->getMaxWindSpeed());
	std::pair<double, double> moistMinMax = std::make_pair(0, planetos->getMaxVapourDensity());
	std::pair<double, double> cWMinMax = std::make_pair(0, planetos->getMaxVapourDensity());
	double h = b.getHeight();
	double dph = PWM::Utils::altitudeAdjustedPressure(h, planetos);
	std::pair<double, double> pMinMax = std::make_pair(90000 - dph, 108000 - dph);

	for (int i = 0; i < b.getObstacles().size(); ++i){
		if (b.getVelocityTheta().getData(i) < vMinMax.first || b.getVelocityTheta().getData(i) > vMinMax.second){
			inRange = false;
			break;
		}
		if (b.getVelocityPhi().getData(i) < vMinMax.first || b.getVelocityPhi().getData(i) > vMinMax.second){
			inRange = false;
			break;
		}
		if (b.getMoisture().getData(i) < moistMinMax.first || b.getMoisture().getData(i) > moistMinMax.second){
			inRange = false;
			break;
		}
		if (b.getCondensedWater().getData(i) < cWMinMax.first || b.getCondensedWater().getData(i) > cWMinMax.second){
			inRange = false;
			break;
		}
		if (b.getPressure().getData(i) < pMinMax.first || b.getPressure().getData(i) > pMinMax.second){
			inRange = false;
			break;
		}
	}
	if (inRange)
		std::cout << "\033[1;32mrandomInit() works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! randomInit() gives wrong results!\033[0m" << std::endl;
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'airLayer'\033[0m:" << std::endl;
	
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
	
	std::cout << "\nTesting of class \033[1;34m'airLayer'\033[0m complete." << std::endl;
}