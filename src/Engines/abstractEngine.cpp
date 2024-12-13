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
 #include "abstractEngine.h"

namespace PWM {
	namespace Engine {
		//Default constructor
        AbstractEngine::AbstractEngine(float t, bool active) : dt(t), execTimePassed(0.0), simTimePassed(0.0), isActive(active){}
		
		//One step through for the engine
		void AbstractEngine::step(){
			if (isActive){
				startComputation();
				step_internal();
				endComputation();
			}
		}
		
		void AbstractEngine::startComputation(){
			execStartPoint = std::chrono::system_clock::now();
		}
		
		void AbstractEngine::endComputation(){
			execTimePassed += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - execStartPoint).count()/1000000000.0;
			simTimePassed += dt;
		}

		float AbstractEngine::getRunTimePassed() const{
			return execTimePassed;
		}

		float AbstractEngine::getSimTimePassed() const{
			return simTimePassed;
		}
		
		const float AbstractEngine::getDt() const{
			return dt;
		}
	}
}
