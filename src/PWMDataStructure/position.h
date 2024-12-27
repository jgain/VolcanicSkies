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
#ifndef PWM_PWMDATASTRUCTURE_POSITION_H
#define PWM_PWMDATASTRUCTURE_POSITION_H

#include <tuple>
namespace PWM{
    namespace PWMDataStructure{
        class Position{
            private:
                double x, y, z;
            public:
                Position(double X = 0, double Y = 0, double Z = 0);
                
                const double getX() const;
                const double getY() const;
                const double getZ() const;
                
                const double getDistanceTo(Position& other) const;
                
                void move(double newX, double newY, double newZ);
                void move(std::tuple<double, double, double>& distance);
                void move(std::tuple<double, double, double>& velocity, double dt);
                void move(Position& newPos);

                bool operator==(Position& other);
                bool operator!=(Position& other);

        };
    }
}
#endif //PWM_PWMDATASTRUCTURE_POSITION_H