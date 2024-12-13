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
#ifndef PWM_PRESSUREWAVEENGINE_H
#define PWM_PRESSUREWAVEENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include "mathUtils.h"
#include "pressureWave.h"
#include <memory>
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename T, typename V> class pressureWaveEngine : public AbstractEngine{
            private:
                //Data structure to hold reference to each of the layers in the simulation
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;

                std::shared_ptr<Model::pressureWave<V>> explosion;
                bool reset, resetOnce;

                float explosionLife = 0.f;

                float pressureWaveLife = 1.5;
            protected:
                void step_internal() override;
            public:
                pressureWaveEngine(float dt, bool active);

                //Function to add an air layer to the engine
                void addLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);
                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;

                //Function to add an explosion to the timeline
                void addExplosion(const std::shared_ptr<Model::pressureWave<V>>& exp);
                float getPressureWaveLife() const;
                void setPressureWaveLife(float newPressureWaveLife);
        };

        template<typename T, typename V>
        inline pressureWaveEngine<T, V>::pressureWaveEngine(float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            resetOnce = reset = false;
        }

        template<typename T, typename V>
        inline void pressureWaveEngine<T, V>::addLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            airLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& pressureWaveEngine<T, V>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename V>
        inline void pressureWaveEngine<T, V>::addExplosion(const std::shared_ptr<Model::pressureWave<V> > &exp){
            explosion = exp;
        }

        template<typename T, typename V>
        inline float pressureWaveEngine<T, V>::getPressureWaveLife() const
        {
            return pressureWaveLife;
        }

        template<typename T, typename V>
        inline void pressureWaveEngine<T, V>::setPressureWaveLife(float newPressureWaveLife)
        {
            pressureWaveLife = newPressureWaveLife;
        }

        template<typename T, typename V>
        inline void pressureWaveEngine<T, V>::step_internal(){
            if (explosion == nullptr)
                return;
            float dT = getDt();
            explosion->updateWave(dT);

            //If explosion isn't started yet, do nothing.
            if (!explosion->isStarted()){
                return;
            }
            this->explosionLife += dT;

            V expMaxLife = explosion->getMaxLife();
            reset = this->explosionLife > (expMaxLife + 10 * dT);

            V zero = 0;
            if (reset){
                for (auto l : airLayers){
                    l->getWavePressure().copy(zero);
                }
                resetOnce = true;
                return;
            }

            if (!resetOnce){
                //Calculate limited area for optimisation
                auto worldLoc = explosion->getCentre();
                V expCentreHeight = explosion->getHeight();
                V radius = explosion->getRadius();
                V expPower = explosion->getInitialPower();
                auto eLoc = std::make_tuple(worldLoc.first, worldLoc.second, expCentreHeight);
                V gridLen = airLayers.at(0)->getWavePressure().gridLength();
                auto offset = airLayers.at(0)->getWavePressure().getOffset();
                auto coordLoc = PWM::Utils::getCoordLoc(worldLoc.first, worldLoc.second, gridLen, offset);
                int coordRadius = std::min((int) (radius / gridLen) + 2, (int) airLayers.at(0)->getWavePressure().getX());

                int lowAirLayerIDX = -1, highAirLayerIDX = 0;
                for (int i = 0; i < airLayers.size(); ++i){
                    if (lowAirLayerIDX < 0)
                        if (airLayers.at(i)->getHeight() > expCentreHeight - radius)
                            lowAirLayerIDX = i;
                    if (airLayers.at(i)->getHeight() < expCentreHeight + radius)
                        highAirLayerIDX = i;
                }

                for (int k = lowAirLayerIDX; k <= highAirLayerIDX; ++k){
                    V airHeight = airLayers.at(k)->getHeight();
                    #pragma omp parallel for
                    for (int i = -coordRadius; i <= coordRadius; ++i){
                        for (int j = -coordRadius; j <= coordRadius; ++j){
//                            if (i == 0 && j == 0 && k == 4)
//                                std::cout << "Debug point reached" << std::endl;
                            auto cellWorldLoc = PWM::Utils::getWorldLoc(coordLoc.first + i, coordLoc.second + j, gridLen, offset);
                            auto cLoc = std::make_tuple(cellWorldLoc.first, cellWorldLoc.second, airHeight);
                            V cellDist = PWM::Utils::calcCartesianDistance<V>(cLoc, eLoc);

                            // Is this cell within the current shockwave radius?
                            if (cellDist < radius){
                                V prevMulVal = airLayers.at(k)->getWavePressure(coordLoc.first + i, coordLoc.second + j);
                                // Check if the radius just reached this cell
                                // If anything in it, the wave has passed and it should now be decreasing.
                                if (prevMulVal != zero){
                                    auto v1 = prevMulVal / timeBitShift;
                                    float v2 = (static_cast<int>(v1)) / 1000.f;
                                    float time = std::abs(v2);
                                    V prevVal = std::abs(std::fmod(prevMulVal, timeBitShift));

                                    /**
                                     * Pressure gradient coefficient dropoff is not linear. Modelled as InitialPower/time.
                                     * Calculate initial gradient coefficient for the cell, then calculate new value.
                                     */
                                    float newTime = time + dT;
                                    V initVal = prevVal * (time + 1);
                                    V newVal = initVal / std::pow((newTime + 1), 1);
                                    V newMulVal = newVal + ((int) (newTime * 1000) * timeBitShift);
                                    if (time > pressureWaveLife)
                                        newMulVal *= -1;
                                    airLayers.at(k)->setWavePressure(coordLoc.first + i, coordLoc.second + j, newMulVal);
                                }
                                // This cell has just entered the shockwave radius
                                else{
                                    float time = 0;
                                    V newVal = expPower / std::pow(cellDist, 1);
                                    V newMulVal = newVal + ((int) (time * 1000) * timeBitShift);
                                    airLayers.at(k)->setWavePressure(coordLoc.first + i, coordLoc.second + j, newMulVal);
                                }
                            }
                        }
                    }
                    #pragma omp barrier
                }
            }
        }
    }
}

#endif // PRESSUREWAVEENGINE_H
