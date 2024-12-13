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
#ifndef PWM_TERRAIN_HEATING_ENGINE_H
#define PWM_TERRAIN_HEATING_ENGINE_H

#include "abstractEngine.h"
#include "atmofuncs.h"
#if __cplusplus > 201703L
    #include <numbers>
#endif
#include "planet.h"
#include "sun.h"
#include "terrain.h"

namespace PWM{
    namespace Engine{
        template<typename T, typename TT, typename V, typename VV> class terrainHeatingEngine : public AbstractEngine{
            private:
                std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> Terrain;
                std::vector<std::shared_ptr<PWM::Model::sun<V>>> sunList;

                V heatLossCoefficient = 0.3;

                void applySunLight(std::shared_ptr<PWM::Model::sun<V>>& s);
            protected:
                void step_internal();
            public:
                terrainHeatingEngine(float dt, bool active);
                terrainHeatingEngine(float dt, bool active, std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t);

                const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& getTerrain() const;
                const std::vector<std::shared_ptr<PWM::Model::sun<V>>>& getSunList() const;
                const V& getHeatLossCoefficient() const;

                void setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t);
                void addSun(std::shared_ptr<PWM::Model::sun<V>>& s);
                void setHeatLossCoefficient(V& hlc);
        };

        template<typename T, typename TT, typename V, typename VV>
        inline terrainHeatingEngine<T, TT, V, VV>::terrainHeatingEngine(float dt, bool active) : AbstractEngine(dt, active){
            Terrain = nullptr;
            sunList = std::vector<std::shared_ptr<PWM::Model::sun<V>>>();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline terrainHeatingEngine<T, TT, V, VV>::terrainHeatingEngine(float dt, bool active, std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t) : AbstractEngine(dt, active), Terrain(t){
            sunList = std::vector<std::shared_ptr<PWM::Model::sun<V>>>();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& terrainHeatingEngine<T, TT, V, VV>::getTerrain() const{
            return Terrain;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<PWM::Model::sun<V>>>& terrainHeatingEngine<T, TT, V, VV>::getSunList() const{
            return sunList;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline  const V& terrainHeatingEngine<T, TT, V, VV>::getHeatLossCoefficient() const{
            return heatLossCoefficient;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrainHeatingEngine<T, TT, V, VV>::setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t){
            Terrain = t;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrainHeatingEngine<T, TT, V, VV>::addSun(std::shared_ptr<PWM::Model::sun<V>>& s){
            sunList.push_back(s);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrainHeatingEngine<T, TT, V, VV>::setHeatLossCoefficient(V& hlc){
            heatLossCoefficient = hlc;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrainHeatingEngine<T, TT, V, VV>::step_internal(){
            if (Terrain == nullptr){
                std::cerr << "Error! Make sure to assign a terrain first! Not stepping through." << std::endl;
                return;
            }
            for (auto s : sunList)
                applySunLight(s);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrainHeatingEngine<T, TT, V, VV>::applySunLight(std::shared_ptr<PWM::Model::sun<V>>& s){
            #if __cplusplus > 201703L
                V pi = std::numbers::pi;
            #else
                V pi = 3.1415926535;
            #endif
                    
            V lumos = s->getPower() * s->getDistance();
            //subsolar point is known
            double sunTheta = PWM::Utils::degToRad<double>(s->getApparentDeclination() + 90);
            double sunPhi = PWM::Utils::degToRad<double>(s->getApparentRightAscension() + 180);
            V atmoCoefficient = 0.7;
            V ratio = Terrain->getPlanet()->getRadius() / Terrain->getPlanet()->getAtmoThickness();
            #pragma omp parallel for
                for (int i = 0; i < Terrain->getTemperature().size(); ++i){
                    //determine sun zenith angle (angular distance from subsolar point)
                    auto coord = Terrain->getTemperature().getCoordinates(i);
                    //Coordinates flipped to 0-pi, 0-2pi instead of -90-90, -180-180 that I use it as elsewhere.
                    //Spherical law of cosines used, requires accuracy of double (not float).
                    double cosOfZ = (std::sin(sunTheta) * std::sin(PWM::Utils::degToRad<double>(coord.getLatitude() + 90))) + (std::cos(sunTheta) * std::cos(PWM::Utils::degToRad<double>(coord.getLatitude()+ 90)) * std::cos(std::abs(sunPhi - PWM::Utils::degToRad<double>(coord.getLongitude() + 180))));
                    
                    V deltaQ = 0;
                    V terTemp = Terrain->getTemperature(i);
                    //Calculate terrain specific heat capacity
                    VV terType = Terrain->getTerrainType(i);
                    V terSpecHeat = Terrain->getPlanet()->getTerrainHeatCapacity().at(terType);
                    if (Terrain->getPlanet()->getTerrainDryType().at(terType))
                        terSpecHeat = (terSpecHeat + (Terrain->getPlanet()->getMoistureHeatCapacity() * Terrain->getMoisture(i))) / (1 + Terrain->getMoisture(i));
                    V deltaTerTemp = 0, newTerTemp = 0;

                    if (!(std::acos(cosOfZ) < (pi / 2))){
                        //Sun is above horizon, do something.
                        //AirMass = sqrt((r cos z)^2 + 2r + 1) - r cos z,
                        //where r = planetRad/planetAtmoThickness, and z is zenith angle
                        //Make assumption of I = 1.1 * I0 (power in vacuum at distance of planet) * 0.7^(AM^0.678)
                        //deltaQ = I * sin(height above horizon) * (1 - albedo)
                        V airMass = std::sqrt(std::pow(ratio * cosOfZ, 2) + (2 * ratio) + 1) - (ratio * cosOfZ);
                        V surfaceLight = 1.1 * lumos * std::pow(atmoCoefficient, std::pow(airMass, 0.678));
                        
                        //Calculate total energy falling on this cell per square meter
                        deltaQ = surfaceLight * cosOfZ * (1 - Terrain->getPlanet()->getTerrainAlbedo().at(terType));

                        //Calculate change in temperature given energy transfer and terrain spec heat
                        deltaTerTemp = (deltaQ * getDt()) / (terSpecHeat * Terrain->getPlanet()->getTerrainDensity().at(terType));
                        
                        //apply change in temp
                        newTerTemp = terTemp + deltaTerTemp;
                        Terrain->setTemperature(i, newTerTemp);
                    }

                    //apply radiation out to space
                    //Grey body radiation with assumption that only 30% of outgoing energy escapes, rest is reflected back or absorbed
                    //An admittedly handwavy way of doing it, but seems to work well enough for now. More tests required.
                    deltaQ = Terrain->getPlanet()->getTerrainEmissivityConstants().at(terType) * stefanBoltzmannConstant * (std::pow(terTemp, 4) - std::pow(2.7, 4)) * 1 * heatLossCoefficient;
                    deltaTerTemp = (deltaQ * getDt()) / (terSpecHeat * Terrain->getPlanet()->getTerrainDensity().at(terType));
                    
                    //apply change in temp
                    newTerTemp = terTemp + deltaTerTemp;
                    Terrain->setTemperature(i, newTerTemp);
                }
            #pragma omp barrier
        }
    }
}
#endif //PWM_TERRAIN_HEATING_ENGINE_H
