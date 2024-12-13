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
#ifndef PWM_ENGINE_ATMOSPHERIC_HEATING_ENGINE_H
#define PWM_ENGINE_ATMOSPHERIC_HEATING_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include "atmofuncs.h"
#include <cmath>
#include <omp.h>
#include "terrain.h"
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename T, typename TT, typename V, typename VV> class atmosphericHeatingEngine : public AbstractEngine{
            private:
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;
                std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> Terrain;

                //The portion of the total atmospheric heating that is a result of surface flux.
                V surfaceFluxRatio;
                
                //The threshold at which radiative exchange takes place
                V threshold;

                void applySurfaceRadiation();
                void sortLayers();
            protected:
                void step_internal();

            public:
                atmosphericHeatingEngine(float dt, bool active, V sfr = 1.5, V thres = 0.0001);

                const V& getSurfaceFluxRatio() const;
                const V& getThreshold() const;
                void setSurfaceFluxRatio(V sfr);
                void setThreshold(V t);

                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& getTerrain() const;
                
                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);
                void setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t);
        };

        template<typename T, typename TT, typename V, typename VV>
        inline atmosphericHeatingEngine<T, TT, V, VV>::atmosphericHeatingEngine(float dt, bool active, V sfr, V thres) : AbstractEngine(dt, active), surfaceFluxRatio(sfr), threshold(thres){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            Terrain = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V& atmosphericHeatingEngine<T, TT, V, VV>::getSurfaceFluxRatio() const{
            return surfaceFluxRatio;
        }
        
        template<typename T, typename TT, typename V, typename VV>
        inline const V& atmosphericHeatingEngine<T, TT, V, VV>::getThreshold() const{
            return threshold;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& atmosphericHeatingEngine<T, TT, V, VV>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::setSurfaceFluxRatio(V sfr){
            surfaceFluxRatio = sfr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::setThreshold(V t){
            threshold = t;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::step_internal(){
            applySurfaceRadiation();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            airLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::sortLayers(){
            std::sort(airLayers.begin(), airLayers.end(), [](std::shared_ptr<PWM::Model::airLayer<T, V>>& l, std::shared_ptr<PWM::Model::airLayer<T, V>>& r){
                return l->getHeight() < r->getHeight();
            });
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& atmosphericHeatingEngine<T, TT, V, VV>::getTerrain() const{
            return Terrain;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t){
            Terrain = t;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void atmosphericHeatingEngine<T, TT, V, VV>::applySurfaceRadiation(){
            V gridlength = Terrain->getElevation().gridLength();
            #pragma omp parallel for
                for (int i = 0; i < Terrain->getTemperature().size(); ++i){
                    V tHeight = Terrain->getElevation(i);
                    int k = 0;
                    for (int j = 0; j < airLayers.size(); ++j){
                        if (airLayers[j]->getHeight() > tHeight){
                            k = j;
                            break;
                        }
                    }
                    auto aL = airLayers[k];
                    
                    //Calculate emmitted radiation, i.e., temp loss
                    V airTemp = aL->getTemperature().getData(i);
                    V terTemp = Terrain->getTemperature(i);
                    auto terType = Terrain->getTerrainType(i);
                    if (std::abs(airTemp - terTemp) < threshold)
                        continue;
                    else if (airTemp < terTemp){
                        V q = Terrain->getPlanet()->getTerrainEmissivityConstants().at(terType) * stefanBoltzmannConstant * (std::pow(terTemp, 4) - std::pow(airTemp, 4)) * 1;
                        V airDensity = PWM::Utils::altitudeAdjustedDensity(aL->getHeight(), Terrain->getPlanet());
                        V terSpecHeat = Terrain->getPlanet()->getTerrainHeatCapacity().at(terType);
                        if (Terrain->getPlanet()->getTerrainDryType().at(terType))
                            terSpecHeat = (terSpecHeat + (aL->getPlanet()->getMoistureHeatCapacity() * Terrain->getMoisture(i))) / (1 + Terrain->getMoisture(i));
                        V deltaTerTemp = (q * getDt()) / (terSpecHeat * Terrain->getPlanet()->getTerrainDensity().at(terType));
                        V newTerTemp = terTemp - deltaTerTemp;
                        Terrain->setTemperature(i, newTerTemp);

                        V deltaAirTemp = ((q * getDt()) / (aL->getPlanet()->getAirHeatCapacity() * airDensity)) * 0.7f; //70% is given as the amount of infrared radiation from the ground that is absorbed by atmosphere. aL->getAbsorptionCoefficient();
                        V newAirTemp = airTemp + deltaAirTemp;
                        aL->setTemperature(i, newAirTemp);
                    }
                }
            #pragma omp barrier
        }
    }
}

#endif //PWM_ENGINE_ATMOSPHERIC_HEATING_ENGINE_H
