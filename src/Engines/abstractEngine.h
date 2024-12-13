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
#ifndef PWM_ABSTRACT_ENGINE_H
#define PWM_ABSTRACT_ENGINE_H

#include <chrono>
namespace PWM{
	namespace Engine{
		class AbstractEngine{
			private:
                //The timestep for this engine in seconds
				float dt;
			
				//Monitoring utilities as used in MWM
				std::chrono::system_clock::time_point execStartPoint;
				float execTimePassed, simTimePassed;
			protected:
				//This function must be overrode by any sub class to provide the actual computation.
                void startComputation();
                void endComputation();
                virtual void step_internal() = 0;

			public:
				//Bool variable to say if this engine is active (e.g., user interaction will start off false until triggered).
				bool isActive;
			
				//Default constructor
				AbstractEngine(float dt = 1.0, bool active = true);

				//Call to perform one step of the engine
				void step();
				//virtual void clear();
				float getRunTimePassed() const;
				float getSimTimePassed() const;
				const float getDt() const;
			/*signals:
				void done();*/
		};
	}
}
#endif //PWM_ABSTRACT_ENGINE_H
