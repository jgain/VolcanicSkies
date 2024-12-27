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
#ifndef PRESSUREWAVE_H
#define PRESSUREWAVE_H

#include "planet.h"
#include <utility>

namespace PWM {
    namespace Model {
        template<typename V> class pressureWave{
        private:
            //The planet this pressure wave/explosion is on
            std::shared_ptr<planet> Planet;
            float startTime; // The time in seconds from start of the simulation to the start of the explosion
            V initialPower;
            V radius;
            std::pair<V, V> centre;
            V height;
            float timeSinceSimStart;

            //The longest that the pressure wave can last
            float maxLife = 50;
            bool started, active;
        public:
            pressureWave();
            pressureWave(std::shared_ptr<planet> P, float startT = 0, V initP = 1000);
            void updateWave(float timeStep);

            float getStartTime() const;
            void setStartTime(float newStartTime);

            V getRadius() const;

            float getTimeSinceStart() const;

            const std::pair<V, V> &getCentre() const;
            const V& getHeight() const;
            void setCentre(const V x, const V y);
            void setCentre(const std::pair<V, V> &newCentre);
            void setHeight(const V h);
            void setPosition(const V x, const V y, const V z);

            V getInitialPower() const;

            bool isActive() const;
            bool isStarted() const;
            float getMaxLife() const;
            void setMaxLife(float newMaxLife);
        };

        template<typename V>
        inline pressureWave<V>::pressureWave()
        {
            Planet = nullptr;
            startTime = -1;
            initialPower = -1;
            radius = 0;
            timeSinceSimStart = 0;
            active = started = false;
        }

        template<typename V>
        inline pressureWave<V>::pressureWave(std::shared_ptr<planet> P, float startT, V initP) : Planet(P), startTime(startT), initialPower(initP){
            radius = 0;
            timeSinceSimStart = 0;
            started = false;
        }

        template<typename V>
        inline void pressureWave<V>::updateWave(float timeStep)
        {
            this->timeSinceSimStart += timeStep;
            if (!started){
                active = started = startTime < timeSinceSimStart;
                return;
            }

            radius += Planet->getSpeedOfSound() * timeStep;

            if (timeSinceSimStart > startTime + maxLife)
                active = false;
        }

        template<typename V>
        inline V pressureWave<V>::getRadius() const
        {
            return radius;
        }

        template<typename V>
        inline float pressureWave<V>::getTimeSinceStart() const
        {
            return timeSinceSimStart;
        }

        template<typename V>
        inline const std::pair<V, V> &pressureWave<V>::getCentre() const
        {
            return centre;
        }

        template<typename V>
        inline const V &pressureWave<V>::getHeight() const
        {
            return this->height;
        }

        template<typename V>
        inline float pressureWave<V>::getMaxLife() const
        {
            return maxLife;
        }

        template<typename V>
        inline void pressureWave<V>::setMaxLife(float newMaxLife)
        {
            maxLife = newMaxLife;
        }

        template<typename V>
        inline bool pressureWave<V>::isStarted() const
        {
            return started;
        }

        template<typename V>
        inline void pressureWave<V>::setCentre(const V x, const V y)
        {
            this->setCentre(std::make_pair(x, y));
        }

        template<typename V>
        inline void pressureWave<V>::setCentre(const std::pair<V, V> &newCentre)
        {
            centre = newCentre;
        }

        template<typename V>
        inline void pressureWave<V>::setHeight(const V h)
        {
            this->height = h;
        }

        template<typename V>
        inline void pressureWave<V>::setPosition(const V x, const V y, const V z)
        {
            setCentre(std::make_pair(x, y));
            setHeight(z);
        }

        template<typename V>
        inline V pressureWave<V>::getInitialPower() const
        {
            return initialPower;
        }

        template<typename V>
        inline bool pressureWave<V>::isActive() const{
            return active;
        }

        template<typename V>
        inline float pressureWave<V>::getStartTime() const
        {
            return startTime;
        }

        template<typename V>
        inline void pressureWave<V>::setStartTime(float newStartTime)
        {
            startTime = newStartTime;
        }
    }
}

#endif // PRESSUREWAVE_H
