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
#ifndef PWM_MOISTURE_UPTAKE_ENGINE_H
#define PWM_MOISTURE_UPTAKE_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include "terrain.h"
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename T, typename TT, typename V, typename VV> class moistureUptakeEngine : public AbstractEngine{
            private:
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;
                std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> Terrain;

                void sortLayers();
                void takeUpMoisture();
            protected:
                void step_internal();
            
            public:
                moistureUptakeEngine(float dt, bool active);

                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& getTerrain() const;
                
                void setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t);
                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l);
        };

        template<typename T, typename TT, typename V, typename VV>
        inline moistureUptakeEngine<T, TT, V, VV>::moistureUptakeEngine(float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            Terrain = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& moistureUptakeEngine<T, TT, V, VV>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& moistureUptakeEngine<T, TT, V, VV>::getTerrain() const{
            return Terrain;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void moistureUptakeEngine<T, TT, V, VV>::setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t){
            Terrain = t;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void moistureUptakeEngine<T, TT, V, VV>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l){
            airLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void moistureUptakeEngine<T, TT, V, VV>::sortLayers(){
            std::sort(airLayers.begin(), airLayers.end(), [](std::shared_ptr<PWM::Model::airLayer<T, V>>& l, std::shared_ptr<PWM::Model::airLayer<T, V>>& r){
                return l->getHeight() < r->getHeight();
            });
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void moistureUptakeEngine<T, TT, V, VV>::step_internal(){
            takeUpMoisture();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void moistureUptakeEngine<T, TT, V, VV>::takeUpMoisture(){
            V latVapHeat = Terrain->getPlanet()->getLatentHeatofVaporisation();
            V moisSpecHeatCap = Terrain->getPlanet()->getMoistureHeatCapacity();
            #pragma omp parallel for
                for (int i = 0; i < Terrain->getElevation().size(); ++i){
                    //figure out which airlayer gets worked on at this point
                    V tHeight = Terrain->getElevation(i);
                    int k = 0;
                    for (int j = 0; j < airLayers.size(); ++j){
                        if (airLayers[j]->getHeight() > tHeight){
                            k = j;
                            break;
                        }
                    }
                    auto aL = airLayers[k];
                    //Do moisture thingies
                    V mois = aL->getMoisture(i);
                    V temp = aL->getTemperature(i);
                    V terMois = Terrain->getMoisture(i);
                    V satVapPres = PWM::Utils::saturationVapourPressure(temp);
                    V atmPres = PWM::Utils::altitudeAdjustedPressure(aL->getHeight(), aL->getPlanet());
                    //V atmPres = l->getPressure().getData(i);
                    V density = PWM::Utils::pressureAdjustedDensity(atmPres, temp, aL->getPlanet());
                    V humidity = mois / density;
                    V partialPressure = atmPres * (humidity / (0.622 + humidity));
                    V satHumidity = (0.622 * satVapPres) / (atmPres - satVapPres);
                    if (atmPres > satVapPres && satVapPres > partialPressure){
                        V diff = satHumidity - humidity;
                        V water = diff * /*transferCoefficient*/ 0.1;
                        if (water > terMois)
                            water = terMois;
                        V newTerMois = terMois - water;
                        V newAirMois = mois + water;
                        aL->setMoisture(i, newAirMois);
                        Terrain->setMoisture(i, newTerMois);
                        V deltaT = water * 0.00001 * latVapHeat / moisSpecHeatCap;
                        V newT = Terrain->getTemperature(i) - deltaT;
                        Terrain->setTemperature(i, newT);
                    }
                    else if (atmPres > partialPressure && partialPressure > satVapPres){
                        V diff = humidity - satHumidity;
                        V water = diff * /*transferCoefficient*/ 0.1;
                        if (water > mois)
                            water = mois;
                        V newTerMois = terMois + water;
                        V newAirMois = mois - water;
                        aL->getMoisture(i, newAirMois);
                        Terrain->setMoisture(i, newTerMois);
                        V deltaT = water * 0.00001 * latVapHeat / moisSpecHeatCap;
                        V newT = Terrain->getTemperature(i) + deltaT;
                        Terrain->setTemperature(i, newT);
                    }
                }
            #pragma omp barrier
        }
    }
}

#endif //PWM_MOISTURE_UPTAKE_ENGINE_H
