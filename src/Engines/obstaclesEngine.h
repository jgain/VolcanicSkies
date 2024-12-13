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
#ifndef PWM_OBSTACLES_ENGINE_H
#define PWM_OBSTACLES_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include "flatStaggeredGrid.h"
#include "terrain.h"
#include <omp.h>
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename T, typename TT, typename V, typename VV> class obstaclesEngine : public AbstractEngine{
            private:
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;

                std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> terrain;

                void sortLayers();
            protected:
                void step_internal();
            public:
                obstaclesEngine(float dt, bool active);

                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l);
                void setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> t);

                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                const std::shared_ptr<PWM::Model::airLayer<T, V>> getAirLayer(const size_t index) const;
                const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> getTerrain() const;
        };

        template<typename T, typename TT, typename V, typename VV>
        inline obstaclesEngine<T, TT, V, VV>::obstaclesEngine(float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void obstaclesEngine<T, TT, V, VV>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l){
            airLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& obstaclesEngine<T, TT, V, VV>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::airLayer<T, V>> obstaclesEngine<T, TT, V, VV>::getAirLayer(const size_t index) const{
            return airLayers.at(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void obstaclesEngine<T, TT, V, VV>::setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> t){
            terrain = t;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> obstaclesEngine<T, TT, V, VV>::getTerrain() const{
            return terrain;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void obstaclesEngine<T, TT, V, VV>::sortLayers(){
            std::sort(airLayers.begin(), airLayers.end(), [](std::shared_ptr<PWM::Model::airLayer<T, V>>& l, std::shared_ptr<PWM::Model::airLayer<T, V>>& r){
                return l->getHeight() < r->getHeight();
            });
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void obstaclesEngine<T, TT, V, VV>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define obstaclesEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void obstaclesEngine<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::step_internal(){
            double zero = 0, one = 1;
            #pragma omp parallel for
                for (int i = 0; i < terrain->getElevation().size(); ++i){
                    for (int j = 0; i < airLayers.size(); ++i){
                        double terrElevation = terrain->getElevation(i);
                        double airHeight = airLayers[i]->getHeight();
                        double airThick = airLayers[i]->getThickness();
                        double sigma = one - std::min(one, std::max(zero, (terrElevation - airHeight) / airThick));
                        double newVelX = airLayers[j]->getVelocityTheta(i) * sigma;
                        double newVelY = airLayers[j]->getVelocityPhi(i) * sigma;
                        airLayers[j]->setVelocityTheta(i, newVelX);
                        airLayers[j]->setVelocityPhi(i, newVelY);
                        if (sigma <= zero){
                            airLayers[j]->setCondensedWater(i, zero);
                            airLayers[j]->setMoisture(i, zero);
                            airLayers[j]->setClouds(i, zero);
                            airLayers[j]->setTemperature(i, terrain->getTemperature(i));
                        }
                    }
                }
            #pragma omp barrier
        }

        template<>
        inline void obstaclesEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::step_internal(){
            double zero = 0, one = 1;
            #pragma omp parallel for
                for (int i = 0; i < terrain->getElevation().size(); ++i){
                    for (int j = 0; i < airLayers.size(); ++i){
                        double terrElevation = terrain->getElevation(i);
                        double airHeight = airLayers[i]->getHeight();
                        double airThick = airLayers[i]->getThickness();
                        double sigmaC = one - std::min(one, std::max(zero, (terrElevation - airHeight) / airThick));
                        if (sigmaC <= zero){
                            airLayers[j]->setCondensedWater(i, zero);
                            airLayers[j]->setMoisture(i, zero);
                            airLayers[j]->setTemperature(i, terrain->getTemperature(i));
                        }

                        auto sigmaThetaLoc = airLayers[j]->getVelocityTheta().getWorldLoc(i);
                        auto sigmaT = one - std::min(one, std::max(zero, (terrain->getElevation().sampleAt(sigmaThetaLoc, airLayers[j]->getObsPtr()) - airHeight) / airThick));
                        double newVelX = airLayers[j]->getVelocityTheta(i) * sigmaT;
                        airLayers[j]->setVelocityTheta(i, newVelX);

                        auto sigmaPhiLoc = airLayers[j]->getVelocityPhi().getWorldLoc(i);
                        auto sigmaP = one - std::min(one, std::max(zero, (terrain->getElevation().sampleAt(sigmaPhiLoc, airLayers[j]->getObsPtr()) - airHeight) / airThick));
                        double newVelY = airLayers[j]->getVelocityPhi(i) * sigmaP;
                        airLayers[j]->setVelocityPhi(i, newVelY);
                    }
                }
            #pragma omp barrier
        }
    }
}

#endif //PWM_OBSTACLES_ENGINE_H
