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
#include <algorithm>
#include <iostream>
#include "planet.h"
#include <string>

using namespace PWM::Model;
void testConstructors(){
    double testRadius = 6371000.f;
    double testRotPeriod = 86400.f;
    double testAveTemp = 288.15f;
    auto a = planet(testRadius, testRotPeriod, testAveTemp);
    bool p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13;
    p1 = a.getRadius() == 6371000.f;
    p2 = a.getRotPeriod() == 86400.f;
    p3 = a.getAveTemp() == 288.15f;

    if (p1 && p2 && p3)
        std::cout << "\033[1;32mDefault constructor works as expected.\033[0m" << std::endl;
    else if (!p1)
        std::cout << "\033[1;31mError! Default constructor does not give expected radius!\033[0m" << std::endl;
    else if (!p2)
        std::cout << "\033[1;31mError! Default constructor does not give expected rotational period!\033[0m" << std::endl;
    else if (!p3)
        std::cout << "\033[1;31mError! Default constructor does not give expected average temperature!\033[0m" << std::endl;
    
    auto b = planet("../resources/PlanetEarth.json");
    p1 = b.getRadius() == 6371000.;
    p2 = b.getRotPeriod() == 86164.1;
    p3 = b.getAveTemp() == 288.15;
    p4 = b.getAverageSealevelPressure() == 101325.;
    p5 = b.getReferenceHeight() == 0.;
    p6 = b.getGravitationalAcceleration() == 9.81;
    p7 = b.getMaxVapourDensity() == 0.294271;
    p8 = b.getMaxWindSpeed() == 40.;
    p9 = b.getMolarAirMass() == 0.0289644;
    p10 = b.getTempLapseRate() == 0.0098;
    p11 = b.getAirHeatCapacity() == 1006;
    p12 = b.getMoistureHeatCapacity() == 4000;
    p13 = b.getLatentHeatofVaporisation() == 2257;
    
    if (p1 && p2 && p3 && p4 && p5 && p6 && p7 && p8 && p9 && p10 && p11 && p12 && p13)
        std::cout << "\033[1;32mDefault JSON constructor works as expected.\033[0m" << std::endl;
    else if (!p1)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected radius!\033[0m" << std::endl;
    else if (!p2)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected rotational period!\033[0m" << std::endl;
    else if (!p3)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected average temperature!\033[0m" << std::endl;
    else if (!p4)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected average sea level pressure!\033[0m" << std::endl;
    else if (!p5)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected reference height!\033[0m" << std::endl;
    else if (!p6)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected gravitational acceleration!\033[0m" << std::endl;
    else if (!p7)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected max vapour density!\033[0m" << std::endl;
    else if (!p8)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected max wind speed!\033[0m" << std::endl;
    else if (!p9)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected molar mass of air!\033[0m" << std::endl;
    else if (!p10)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected temperature lapse rate!\033[0m" << std::endl;
    else if (!p11)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected air specific heat capacity!\033[0m" << std::endl;
    else if (!p12)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected moisture specific heat capacity!\033[0m" << std::endl;
    else if (!p13)
        std::cout << "\033[1;31mError! Default JSON constructor does not give expected latent heat of vaporisation!\033[0m" << std::endl;
    
}

void testComparators(){
	auto a = planet("../resources/PlanetEarth.json");
    auto b = planet("../resources/PlanetEarth.json");
    if (a == b)
        std::cout << "\033[1;32mEquality operator works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator does not perform as expected!\033[0m" << std::endl;

    double r = 10;
    b.setRadius(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing radii.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing radii!\033[0m" << std::endl;

    b = a;
    r = 1440;
    b.setRotPeriod(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing rotational periods.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing rotational periods!\033[0m" << std::endl;

    b = a;
    r = 300;
    b.setAveTemp(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing average temps.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing average temps!\033[0m" << std::endl;

    b = a;
    r = 98000;
    b.setAverageSealevelPressure(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing average sea level pressures.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing average sea level pressures!\033[0m" << std::endl;

    b = a;
    r = 9;
    b.setReferenceHeight(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing reference heights.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing reference heights!\033[0m" << std::endl;

    b = a;
    r = 6.5;
    b.setGravitationalAcceleration(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing gravitational accelerations.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing gravitational accelerations!\033[0m" << std::endl;

    b = a;
    r = 0.6;
    b.setMaxVapourDensity(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing max vapour densities.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing max vapour densities!\033[0m" << std::endl;

    b = a;
    r = 90;
    b.setMaxWindSpeed(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing max wind speeds.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing max wind speeds!\033[0m" << std::endl;

    b = a;
    r = 0.3;
    b.setMolarAirMass(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing molar air masses.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing molar air masses!\033[0m" << std::endl;

    b = a;
    r = 5.5;
    b.setTempLapseRate(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing temperature lapse rates.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing temperature lapse rates!\033[0m" << std::endl;
    
    b = a;
    auto x = b.getTerrainTypes();
    x.push_back("FalseTerrain");
    b.setTerrainTypes(x);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing terrain type vectors.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing terrain type vectors!\033[0m" << std::endl;

    b = a;
    auto t = b.getTerrainDensity();
    t.at("WATER") = 984;
    b.setTerrainDensity(t);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing terrain density maps.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing terrain density maps!\033[0m" << std::endl;

    b = a;
    b.setAirHeatCap(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing air heat capacities.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing air heat capacities!\033[0m" << std::endl;

    b = a;
    b.setMoistureHeatCap(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing moisture heat capacities.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing moisture heat capacities!\033[0m" << std::endl;

    b = a;
    t = b.getTerrainHeatCapacity();
    t.at("WATER") = 48;
    b.setTerrainHeatCapacity(t);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing terrain heat capacity maps.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing terrain heat capacity maps!\033[0m" << std::endl;

    b = a;
    t = b.getTerrainAlbedo();
    t.at("WATER") = 0.2;
    b.setTerrainAlbedo(t);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing terrain albedo maps.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing terrain albedo maps!\033[0m" << std::endl;

    b = a;
    t = b.getTerrainEmissivityConstants();
    t.at("ICE") = 0.8;
    b.setTerrainEmissivityConstants(t);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing terrain emissivity maps.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing terrain emissivity maps!\033[0m" << std::endl;
    
    b = a;
    auto t2 = b.getTerrainDryType();
    t2.at("ICE") = true;
    b.setTerrainDryType(t2);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing terrain dry maps.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing terrain dryy maps!\033[0m" << std::endl;
    
    b = a;
    b.setLatentHeatofVaporisation(r);
    if (a != b)
        std::cout << "\033[1;32mEquality operator correctly rejects differing latent heats of vaporisation.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Equality operator accepts differing latent heats of vaporisation!\033[0m" << std::endl;

}

void testGetters(){
	auto a = planet("../resources/PlanetEarth.json");
	
    if (a.getRadius() == 6371000.)
        std::cout << "\033[1;32mgetRadius() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getRadius() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getRotPeriod() == 86164.1)
        std::cout << "\033[1;32mgetRotPeriod() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getRotPeriod() does not perform as expected!\033[0m" << std::endl;

    if (a.getAveTemp() == 288.15)
        std::cout << "\033[1;32mgetAveTemp() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAveTemp() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getAverageSealevelPressure() == 101325.)
        std::cout << "\033[1;32mgetAverageSealevelPressure() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAverageSealevelPressure() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getReferenceHeight() == 0.)
        std::cout << "\033[1;32mgetReferenceHeight() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getReferenceHeight() does not perform as expected!\033[0m" << std::endl;

    if (a.getGravitationalAcceleration() == 9.81)
        std::cout << "\033[1;32mgetGravitationalAcceleration() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getGravitationalAcceleration() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getMaxVapourDensity() == 0.294271)
        std::cout << "\033[1;32mgetMaxVapourDensity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getMaxVapourDensity() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getMaxWindSpeed() == 40.)
        std::cout << "\033[1;32mgetMaxWindSpeed() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getMaxWindSpeed() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getMolarAirMass() == 0.0289644)
        std::cout << "\033[1;32mgetMolarAirMass() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getMolarAirMass() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getTempLapseRate() == 0.0098)
        std::cout << "\033[1;32mgetTempLapseRate() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTempLapseRate() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getAirHeatCapacity() == 1006)
        std::cout << "\033[1;32mgetAirHeatCapacity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAirHeatCapacity() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getMoistureHeatCapacity() == 4000)
        std::cout << "\033[1;32mgetMoistureHeatCapacity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getMoistureHeatCapacity() does not perform as expected!\033[0m" << std::endl;
    
    auto x = a.getTerrainTypes();
    if (std::find(x.begin(), x.end(), "ICE") != x.end())
        std::cout << "\033[1;32mgetTerrainTypes() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrainTypes() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getTerrainDensity().at("WATER") == 1000)
        std::cout << "\033[1;32mgetTerrainDensity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrainDensity() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getTerrainHeatCapacity().at("ICE") == 2104)
        std::cout << "\033[1;32mgetTerrainHeatCapacity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrainHeatCapacity() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getTerrainAlbedo().at("ICE") == 0.4)
        std::cout << "\033[1;32mgetTerrainAlbedo() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrainAlbedo() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getTerrainEmissivityConstants().at("ICE") == 0.98)
        std::cout << "\033[1;32mgetTerrainEmissivityConstants() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrainEmissivityConstants() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getTerrainDryType().at("ICE") == false)
        std::cout << "\033[1;32mgetTerrainDryType() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getTerrainDryType() does not perform as expected!\033[0m" << std::endl;
    
    if (a.getLatentHeatofVaporisation() == 2257)
        std::cout << "\033[1;32mgetLatentHeatofVaporisation() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getLatentHeatofVaporisation() does not perform as expected!\033[0m" << std::endl;
}

void testSetters(){
	auto a = planet("../resources/PlanetEarth.json");
	
    double r1 = 5.;
    a.setRadius(r1);
    if (a.getRadius() == 5.)
        std::cout << "\033[32msetRadius() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setRadius() does not perform as expected!\033[0m" << std::endl;
    
    a.setRotPeriod(r1);
    if (a.getRotPeriod() == 5.)
        std::cout << "\033[32msetRotPeriod() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setRotPeriod() does not perform as expected!\033[0m" << std::endl;

    a.setAveTemp(r1);
    if (a.getAveTemp() == 5.)
        std::cout << "\033[32msetAveTemp() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setAveTemp() does not perform as expected!\033[0m" << std::endl;

    a.setAverageSealevelPressure(r1);
    if (a.getAverageSealevelPressure() == 5.)
        std::cout << "\033[32msetAverageSealevelPressure() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setAverageSealevelPressure() does not perform as expected!\033[0m" << std::endl;

    a.setReferenceHeight(r1);
    if (a.getReferenceHeight() == 5.)
        std::cout << "\033[32msetReferenceHeight() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setReferenceHeight() does not perform as expected!\033[0m" << std::endl;

    a.setGravitationalAcceleration(r1);
    if (a.getGravitationalAcceleration() == 5.)
        std::cout << "\033[32msetGravitationalAcceleration() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setGravitationalAcceleration() does not perform as expected!\033[0m" << std::endl;

    a.setMaxVapourDensity(r1);
    if (a.getMaxVapourDensity() == 5.)
        std::cout << "\033[32msetMaxVapourDensity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setMaxVapourDensity() does not perform as expected!\033[0m" << std::endl;

    a.setMaxWindSpeed(r1);
    if (a.getMaxWindSpeed() == 5.)
        std::cout << "\033[32msetMaxWindSpeed() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setMaxWindSpeed() does not perform as expected!\033[0m" << std::endl;

    a.setMolarAirMass(r1);
    if (a.getMolarAirMass() == 5.)
        std::cout << "\033[32msetMolarAirMass() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setMolarAirMass() does not perform as expected!\033[0m" << std::endl;

    a.setTempLapseRate(r1);
    if (a.getTempLapseRate() == 5.)
        std::cout << "\033[32msetTempLapseRate() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[31mError! setTempLapseRate() does not perform as expected!\033[0m" << std::endl;
    
    a.setAirHeatCap(r1);
    if (a.getAirHeatCapacity() == 5.)
        std::cout << "\033[32msetAirHeatCapacity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setAirHeatCapacity() does not perform as expected!\033[0m" << std::endl;
    
    a.setMoistureHeatCap(r1);
    if (a.getMoistureHeatCapacity() == 5.)
        std::cout << "\033[32msetMoistureHeatCapacity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setMoistureHeatCapacity() does not perform as expected!\033[0m" << std::endl;
    
    auto x(a.getTerrainTypes());
    std::string t1 = "FAKE";
    x.push_back(t1);
    a.setTerrainTypes(x);
    if (std::find(a.getTerrainTypes().begin(), a.getTerrainTypes().end(), t1) != a.getTerrainTypes().end())
        std::cout << "\033[32msetTerrainTypes() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrainTypes() does not perform as expected!\033[0m" << std::endl;
    
    auto t = a.getTerrainDensity();
    t.insert({t1, 1});
    a.setTerrainDensity(t);
    if (a.getTerrainDensity().at(t1) == 1)
        std::cout << "\033[32msetTerrainDensity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrainDensity() does not perform as expected!\033[0m" << std::endl;
    
    t = a.getTerrainHeatCapacity();
    t.insert({t1, 1});
    a.setTerrainHeatCapacity(t);
    if (a.getTerrainHeatCapacity().at(t1) == 1)
        std::cout << "\033[32msetTerrainHeatCapacity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrainHeatCapacity() does not perform as expected!\033[0m" << std::endl;
    
    t = a.getTerrainAlbedo();
    t.insert({t1, 1});
    a.setTerrainAlbedo(t);
    if (a.getTerrainAlbedo().at(t1) == 1)
        std::cout << "\033[32msetTerrainAlbedo() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrainAlbedo() does not perform as expected!\033[0m" << std::endl;
    
    t = a.getTerrainEmissivityConstants();
    t.insert({t1, 1});
    a.setTerrainEmissivityConstants(t);
    if (a.getTerrainEmissivityConstants().at(t1) == 1)
        std::cout << "\033[32msetTerrainEmissivityConstants() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrainEmissivityConstants() does not perform as expected!\033[0m" << std::endl;
    
    auto t5 = a.getTerrainDryType();
    t5.insert({t1, false});
    a.setTerrainDryType(t5);
    if (a.getTerrainDryType().at(t1) == false)
        std::cout << "\033[32msetTerrainDryType() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setTerrainDryType() does not perform as expected!\033[0m" << std::endl;
    
    a.setLatentHeatofVaporisation(r1);
    if (a.getLatentHeatofVaporisation() == 5.)
        std::cout << "\033[32msetLatentHeatofVaporisation() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! setLatentHeatofVaporisation() does not perform as expected!\033[0m" << std::endl;
    
}

void testMiscellany(){
	auto a = planet("../resources/PlanetEarth.json");
    double b = (2 * 3.1415926535) / 86164.1;
    if (a.getAngularVelocity() == b)
        std::cout << "\033[1;32mgetAngularVelocity() works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! getAngularVelocity() does not perform as expected!\033[0m" << std::endl;
    
    std::cout << "The following two tests should fail safely, catching exceptions thrown by the JSON library." << std::endl;
    std::cout << "Testing empty string..." << std::endl;
    auto c = planet("");

    std::cout << "Testing wrongly formatted/incomplete JSON file..." << std::endl;
    auto d = planet("../resources/WrongFormattedPlanet.json");
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'Planet'\033[0m:" << std::endl;
	
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
	
	std::cout << "\nTesting of class \033[1;34m'Planet'\033[0m complete." << std::endl;
}