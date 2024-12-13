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
// Wrapper class to deal with all coordinate data (latitude/longitude)
// Cilliers Pretorius
// 07 May 2021

#ifndef PWM_PWMDATASTRUCTURE_COORDINATE_H
#define PWM_PWMDATASTRUCTURE_COORDINATE_H
namespace PWM {
	namespace PWMDataStructure{
		class Coordinate{
			private:
				double latitude, longitude;
				static constexpr double approxLimit = 0.001;
			public:
				void setLatitude(const double lat);
				void setLongitude(const double longi);
				double getLatitude() const;
				double getLongitude() const;
				Coordinate(double lat, double longi);
				/*Coordinate(const Coordinate& other);
				Coordinate & operator=(const Coordinate& other);*/
				bool operator==(const Coordinate& other) const;
				bool operator!=(const Coordinate& other) const;
				bool approxEquals(const Coordinate& other) const;
		};
	}
}

#endif //PWM_PWMDATASTRUCTURE_COORDINATE_H