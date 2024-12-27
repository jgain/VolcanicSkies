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
#include "Coordinate.h"
#include <cmath>

namespace PWM{
	namespace PWMDataStructure{
		Coordinate::Coordinate(double lat, double longi): latitude(lat), longitude(longi) {
		}
		
		/*Coordinate::Coordinate(const Coordinate& other): latitude(other.getLatitude()), longitude(other.getLongitude()) {
		}
		
		Coordinate& Coordinate::operator= (const Coordinate& other){
			setLatitude(other.getLatitude());
			setLongitude(other.getLongitude());
			return *this;
		}*/
		
		bool Coordinate::operator==(const Coordinate& other) const{
			return (getLatitude() == other.getLatitude()) && (getLongitude() == other.getLongitude());
		}

		bool Coordinate::operator!=(const Coordinate& other) const{
			return !(*this == other);
		}

		bool Coordinate::approxEquals(const Coordinate& other) const{
			return (std::abs(getLatitude() - other.getLatitude()) < approxLimit) && (std::abs(getLongitude() - other.getLongitude()) < approxLimit);
		}
		
		double Coordinate::getLatitude() const{
			return this->latitude;
		}
		
		double Coordinate::getLongitude() const{
			return this->longitude;
		}
		
		void Coordinate::setLatitude(const double lat){
			this->latitude = lat;
		}
		void Coordinate::setLongitude(const double longi){
			this->longitude = longi;
		}
				
	}
}