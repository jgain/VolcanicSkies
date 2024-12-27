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
#ifndef PWM_PRECIPITATION_ENGINE_H
#define PWM_PRECIPITATION_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include "convectionLayer.h"
#include "terrain.h"

namespace PWM {
    namespace Engine {
        template <typename T, typename TT, typename V, typename VV> class precipitationEngine : public AbstractEngine{
            private:
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;
                std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>> convectionLayers;
                std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>> terr;

                void sortLayers();

                // Factor applied to vertical velocity for coalescence factor calculation
                V vertVelCoalFactor = 0.1;

                // Factor applied to horizontal velocity for coalescence factor calculation
                V horVelCoalFactor = 0.5;

                // Factor applied to calculate rainfall potential
                V rainLiftFactor = 0.5;

                // Factor applied to calculate ashfall potential
                V ashLiftFactor = 0.001;

                // Factor applied to simulate electrostatic attraction and surface tension between water molecules
                V baseRainFallResist = 2;

                // Factor applied to simulate electrostatic attraction between ash particles
                V baseAshFallResist = 0.02;

                // Factor applied to rain falling on ground to determine absorption
                V groundAbsorbFactor = 1;

                // Factor applied to raindrops falling out of a vertical layer, roughly approximating terminal velocity
                V rainVertFallFactor = 10;

                // Factor applied to ash falling out of a vertical layer, roughly approximating terminal velocity
                V ashVertFallFactor = 100;

                void rainFall();

                void ashFall();
            protected:
                void step_internal() override;
            public:
                precipitationEngine(float dt, bool active);
                precipitationEngine(float dt, bool active, V rLiftFac = 0.05, V horVelCoFac = 1, V vertVelCoFac = 1);
                precipitationEngine(float dt, bool active, V rLiftFac, V aLiftFac, V horVelCoFac, V vertVelCoFac);

                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& getConvectionLayers() const;
                const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& getTerrain() const;

                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);
                void addConvectionLayer(std::shared_ptr<PWM::Model::convectionLayer<T, V>>& l);
                void setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t);
        };

        template<typename T, typename TT, typename V, typename VV>
        inline precipitationEngine<T, TT, V, VV>::precipitationEngine(float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            terr = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline precipitationEngine<T, TT, V, VV>::precipitationEngine(float dt, bool active, V rLiftFac, V horVelCoFac, V vertVelCoFac) : AbstractEngine(dt, active), rainLiftFactor(rLiftFac), horVelCoalFactor(horVelCoFac), vertVelCoalFactor(vertVelCoFac){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            terr = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline precipitationEngine<T, TT, V, VV>::precipitationEngine(float dt, bool active, V rLiftFac, V aLiftFac, V horVelCoFac, V vertVelCoFac) : AbstractEngine(dt, active), rainLiftFactor(rLiftFac), ashLiftFactor(aLiftFac), horVelCoalFactor(horVelCoFac), vertVelCoalFactor(vertVelCoFac){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            terr = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& precipitationEngine<T, TT, V, VV>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& precipitationEngine<T, TT, V, VV>::getConvectionLayers() const{
            return convectionLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& precipitationEngine<T, TT, V, VV>::getTerrain() const{
            return terr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            airLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::addConvectionLayer(std::shared_ptr<PWM::Model::convectionLayer<T, V>>& l){
            convectionLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::setTerrain(std::shared_ptr<PWM::Model::terrain<T, TT, V, VV>>& t){
            terr = t;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::sortLayers(){
            std::sort(airLayers.begin(), airLayers.end());
            std::sort(convectionLayers.begin(), convectionLayers.end());
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::rainFall(){
            float deltaT = this->getDt();
            for (int i = 0; i < convectionLayers.size(); ++i){
                std::shared_ptr<PWM::Model::airLayer<T, V>> botLayer = airLayers[i];
                std::shared_ptr<PWM::Model::airLayer<T, V>> topLayer = airLayers[i + 1];
                std::shared_ptr<PWM::Model::convectionLayer<T, V>> midLayer = convectionLayers[i];
                V heightFactor = 1 / midLayer->getThickness();
                int size = botLayer->getObstacles().size();
                V coalescenceAvg = 0.;
                #pragma omp parallel for
                for (int j = 0; j < size; ++j){
                    if (topLayer->getObstacles(j))
                        continue;

                    V vertVel = midLayer->getVerticalVelocity(j);

                    V velX = topLayer->getVelocityPhi(j);
                    V velY = topLayer->getVelocityTheta(j);
                    V horVel = std::sqrt(std::pow(velX, 2) + std::pow(velY, 2));
                    V coalescenceFactor = std::abs(vertVel) * vertVelCoalFactor + horVel * horVelCoalFactor;
                    if(coalescenceFactor > 1.0)
                        coalescenceFactor = 1.0;
                    coalescenceAvg += coalescenceFactor;

                    V botWater = botLayer->getCondensedWater(j);
                    V topWater = topLayer->getCondensedWater(j);

                    V currRF = midLayer->getRainfall(j);
                    V rfOut = currRF * rainVertFallFactor * deltaT * heightFactor;
                    if (rfOut > currRF)
                        rfOut = currRF;

                    //First, check if terrain is in bottom layer. If so, sploosh.
                    if (terr->getElevation(j) >= botLayer->getLayerBot()){
                        V terrMoist = terr->getMoisture(j);
                        V newTM = terrMoist + rfOut;
                        terr->setMoisture(j, newTM);
                    }
                    //Else, just drop the water into the layer below (possibly through? test).
                    else{
                        V newBotWater = botWater + rfOut;
                        botLayer->setCondensedWater(j, newBotWater);
                    }

                    
                    // V potentialRF = topWater * coalescenceFactor * deltaT * heightFactor;
                    // JG - test of rain initiation should not depend on deltT and heightFactor
                    // Also coalescenceFactor needs to be capped at 1.0
                    
                    V potentialRF = topWater * coalescenceFactor;
                    V fallResist = baseRainFallResist + vertVel * rainLiftFactor;
                    bool fall = fallResist < potentialRF;
                    V rfFall = 0.;
                    if (fall){
                        rfFall = potentialRF * rainVertFallFactor * deltaT * heightFactor;
                        if (rfFall > topWater)
                            rfFall = topWater;
                        V newTopWater = topWater - rfFall;
                        topLayer->setCondensedWater(j, newTopWater);
                    }

                    V newRF = currRF - rfOut + rfFall;
                    midLayer->setRainfall(j, newRF);
                }
                 #pragma omp barrier
                coalescenceAvg /= (float) size;
                // std::cerr << "Average coallescence for layer " << i << " = " << coalescenceAvg << std::endl;
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::ashFall(){
            float deltaT = this->getDt();
            for (int i = 0; i < convectionLayers.size(); ++i){
                std::shared_ptr<PWM::Model::airLayer<T, V>> botLayer = airLayers[i];
                std::shared_ptr<PWM::Model::airLayer<T, V>> topLayer = airLayers[i + 1];
                std::shared_ptr<PWM::Model::convectionLayer<T, V>> midLayer = convectionLayers[i];
                V heightFactor = 1 / midLayer->getThickness();
                int size = botLayer->getObstacles().size();
                #pragma omp parallel for
                for (int j = 0; j < size; ++j){
                    if (topLayer->getObstacles(j))
                        continue;

                    V vertVel = midLayer->getVerticalVelocity(j);

                    V velX = topLayer->getVelocityPhi(j);
                    V velY = topLayer->getVelocityTheta(j);
                    V horVel = std::sqrt(std::pow(velX, 2) + std::pow(velY, 2));
                    V coalescenceFactor = std::abs(vertVel) * vertVelCoalFactor + horVel * horVelCoalFactor;
                    if(coalescenceFactor > 1.0)
                        coalescenceFactor = 1.0;

                    V botParts = botLayer->getParticulates(j);
                    V topParts = topLayer->getParticulates(j);

                    V currRF = midLayer->getAshfall(j);
                    V rfOut = currRF * ashVertFallFactor * (deltaT * heightFactor);
                    if (rfOut > currRF)
                        rfOut = currRF;

                    //First, check if terrain is in bottom layer. If so, drop the ash on the ground.
                    if (!(terr->getElevation(j) >= botLayer->getLayerBot())){
                        terr->addAsh(j, rfOut);
                    }
                    //Else, just drop the water into the layer below (possibly through? test).
                    else{
                        V newBotParts = botParts + rfOut;
                        botLayer->setParticulates(j, newBotParts);
                    }

                    // V potentialAF = topParts * coalescenceFactor * deltaT * heightFactor;
                    // JG - test of rain initiation should not depend on deltT and heightFactor
                    // Also coalescenceFactor needs to be capped at 1.0
                    
                    V potentialAF = topParts * coalescenceFactor;

                    V fallResist = baseAshFallResist + vertVel * ashLiftFactor;
                    bool fall = fallResist < potentialAF;
                    V afFall = 0.;
                    if (fall){
                        afFall = potentialAF * ashVertFallFactor * deltaT * heightFactor;
                        if (afFall > topParts)
                            afFall = topParts;
                        V newTopParts = topParts - afFall;
                        topLayer->setParticulates(j, newTopParts);
                    }

                    V newRF = currRF - rfOut + afFall;
                    midLayer->setAshfall(j, newRF);
                }
                #pragma omp barrier
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void precipitationEngine<T, TT, V, VV>::step_internal(){
            rainFall();
            ashFall();
        }
    }
}

#endif // PWM_PRECIPITATION_ENGINE_H
