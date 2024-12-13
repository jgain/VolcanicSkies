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
#ifndef PWM_ENGINE_VERTICAL_COUPLING_ENGINE_H
#define PWM_ENGINE_VERTICAL_COUPLING_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include <algorithm>
#include "atmofuncs.h"
#include "convectionLayer.h"
#include "flatStaggeredGrid.h"
#include "mathUtils.h"
#include "square2DArray.h"
#include <omp.h>
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename T, typename V> class verticalCouplingEngine : public AbstractEngine{
            private:
                void calcConvectiveUplift(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer);
                void calcPressureUplift(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer);
                void transportFromUplift(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer);
                void upliftInOutflowApproximation(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer);
                void upliftInOutflowApproximation();
                void smoothVertVelocities();

                V convectCouplingCoeff = 0.001;
                V pressureCouplingCoeff = 0.05;
                V pressureUpliftCoeff = 0.005;
                V upliftTransportCoeff = 10.;
                V moveAverageWindowRadius = 3;
                bool performPressureUplift;
                bool performPressureSmooth = false;
                bool performVertVelSmooth = true;

                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;
                std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>> convectionLayers;
                void sortLayers();

                std::shared_ptr<T> bufferCloudB;
                std::shared_ptr<T> bufferCloudT;
                std::shared_ptr<T> bufferMoistureB;
                std::shared_ptr<T> bufferMoistureT;
                std::shared_ptr<T> bufferPressureB;
                std::shared_ptr<T> bufferPressureT;
                std::shared_ptr<T> bufferTemperatureB;
                std::shared_ptr<T> bufferTemperatureT;
//                std::shared_ptr<T> bufferVelThetaB;
//                std::shared_ptr<T> bufferVelThetaT;
//                std::shared_ptr<T> bufferVelPhiB;
                //                std::shared_ptr<T> bufferVelPhiT;
            protected:
                void step_internal() override;

            public:
                verticalCouplingEngine(float dt, bool active, bool presUp = false);
                verticalCouplingEngine(size_t theta, size_t phi, float dt, bool active, bool presUp = false);
                verticalCouplingEngine(const size_t faceWidth, float dt, bool active, bool presUp = false);
                verticalCouplingEngine(const size_t xWidth, const size_t yHeight, const V xSize, const V ySize, float dt, bool active, bool presUp = false);

                const V& getConvectCouplingCoeff() const;
                const V& getPressureCouplingCoeff() const;
                const V& getPressureUpliftCoeff() const;

                const bool getPressureUplift() const;
                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& getConvectionLayers() const;

                void setConvectCouplingCoeff(const V& v);
                void setPressureCouplingCoeff(const V& v);
                void setPressureUpliftCoeff(const V& v);

                void setPressureUplift(const bool v);
                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);
                void addConvectionLayer(std::shared_ptr<PWM::Model::convectionLayer<T, V>>& l);

                V debugConvectiveMax = 120;
        };

        template<typename T, typename V>
        inline verticalCouplingEngine<T, V>::verticalCouplingEngine(float dt, bool active, bool presUp) : AbstractEngine(dt, active), performPressureUplift(presUp) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            bufferCloudB = std::make_shared<T>();
            bufferCloudT = std::make_shared<T>();
            bufferMoistureB = std::make_shared<T>();
            bufferMoistureT = std::make_shared<T>();
            bufferPressureB = std::make_shared<T>();
            bufferPressureT = std::make_shared<T>();
            bufferTemperatureB = std::make_shared<T>();
            bufferTemperatureT = std::make_shared<T>();
//            bufferVelThetaB = std::make_shared<T>();
//            bufferVelThetaT = std::make_shared<T>();
//            bufferVelPhiB = std::make_shared<T>();
//            bufferVelPhiT = std::make_shared<T>();
        }

        template<typename T, typename V>
        inline verticalCouplingEngine<T, V>::verticalCouplingEngine(size_t theta, size_t phi, float dt, bool active, bool presUp) : AbstractEngine(dt, active), performPressureUplift(presUp) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            bufferCloudB = std::make_shared<T>(theta, phi);
            bufferCloudT = std::make_shared<T>(theta, phi);
            bufferMoistureB = std::make_shared<T>(theta, phi);
            bufferMoistureT = std::make_shared<T>(theta, phi);
            bufferPressureB = std::make_shared<T>(theta, phi);
            bufferPressureT = std::make_shared<T>(theta, phi);
            bufferTemperatureB = std::make_shared<T>(theta, phi);
            bufferTemperatureT = std::make_shared<T>(theta, phi);
//            bufferVelThetaB = std::make_shared<T>(theta + 1, phi + 1);
//            bufferVelThetaT = std::make_shared<T>(theta + 1, phi + 1);
//            bufferVelPhiB = std::make_shared<T>(theta + 1, phi + 1);
//            bufferVelPhiT = std::make_shared<T>(theta + 1, phi + 1);
        }

        template<typename T, typename V>
        inline verticalCouplingEngine<T, V>::verticalCouplingEngine(const size_t faceWidth, float dt, bool active, bool presUp) : AbstractEngine(dt, active), performPressureUplift(presUp) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            bufferCloudB = std::make_shared<T>(faceWidth);
            bufferCloudT = std::make_shared<T>(faceWidth);
            bufferMoistureB = std::make_shared<T>(faceWidth);
            bufferMoistureT = std::make_shared<T>(faceWidth);
            bufferPressureB = std::make_shared<T>(faceWidth);
            bufferPressureT = std::make_shared<T>(faceWidth);
            bufferTemperatureB = std::make_shared<T>(faceWidth);
            bufferTemperatureT = std::make_shared<T>(faceWidth);
//            bufferVelThetaB = std::make_shared<T>(faceWidth);
//            bufferVelThetaT = std::make_shared<T>(faceWidth);
//            bufferVelPhiB = std::make_shared<T>(faceWidth);
//            bufferVelPhiT = std::make_shared<T>(faceWidth);
        }

        template<typename T, typename V>
        inline verticalCouplingEngine<T, V>::verticalCouplingEngine(const size_t xWidth, const size_t yHeight, const V xSize, const V ySize, float dt, bool active, bool presUp) : AbstractEngine(dt, active), performPressureUplift(presUp) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>();
            bufferCloudB = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferCloudT = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferMoistureB = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferMoistureT = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferPressureB = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferPressureT = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferTemperatureB = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            bufferTemperatureT = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);

            int moveWindowAve = std::round(xWidth / 150.);
            this->moveAverageWindowRadius = moveWindowAve;
        }

        template<typename T, typename V>
        inline const V& verticalCouplingEngine<T, V>::getConvectCouplingCoeff() const{
            return convectCouplingCoeff;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::setConvectCouplingCoeff(const V& v){
            convectCouplingCoeff = v;
        }

        template<typename T, typename V>
        inline const V& verticalCouplingEngine<T, V>::getPressureUpliftCoeff() const{
            return pressureUpliftCoeff;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::setPressureUpliftCoeff(const V& v){
            pressureUpliftCoeff = v;
        }

        template<typename T, typename V>
        inline const V& verticalCouplingEngine<T, V>::getPressureCouplingCoeff() const{
            return pressureCouplingCoeff;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::setPressureCouplingCoeff(const V& v){
            pressureCouplingCoeff = v;
        }

        template<typename T, typename V>
        inline const bool verticalCouplingEngine<T, V>::getPressureUplift() const{
            return performPressureUplift;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::setPressureUplift(const bool v){
            performPressureUplift = v;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            airLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::addConvectionLayer(std::shared_ptr<PWM::Model::convectionLayer<T, V>>& l){
            convectionLayers.push_back(l);
            sortLayers();
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& verticalCouplingEngine<T, V>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& verticalCouplingEngine<T, V>::getConvectionLayers() const{
            return convectionLayers;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::sortLayers(){
            std::sort(airLayers.begin(), airLayers.end());
            std::sort(convectionLayers.begin(), convectionLayers.end());
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::calcConvectiveUplift(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer){
            if (bottomLayer->getObstacles().size() != topLayer->getObstacles().size() || bottomLayer->getObstacles().size() != middleLayer->getVerticalVelocities().size()){
                std::cerr << "Mismatch in grid sizes at verticalCouplingEngine<T, V>::convect()!" << std::endl;
                return;
            }

            V vSolid = 0.0;//Since mountains don't walk or jump, always 0.
            V Cp = bottomLayer->getPlanet()->getAirHeatCapacity();
            V g = bottomLayer->getPlanet()->getGravitationalAcceleration();
            V latVap = bottomLayer->getPlanet()->getLatentHeatofVaporisation();

            //Altitudes
            V hBot = bottomLayer->getHeight();
            V hTop = topLayer->getHeight();

            //Air pressures as calculated by barometric formula
            //Why not use the pressure from airLayer model itself? Edit much later: Because it is stupidly small values, not suitable for MSE.
            V pBot = PWM::Utils::altitudeAdjustedPressure(bottomLayer->getHeight(), bottomLayer->getPlanet());
            V pTop = PWM::Utils::altitudeAdjustedPressure(topLayer->getHeight(), topLayer->getPlanet());

            V dryLapseRate = g / Cp;
            V wetLapseRate = 5. / 1000.; //Approximation, actual value depends on temperature and pressure

            V heightScaleFactor = (1.0 / bottomLayer->getObstacles().gridLength()) / (hTop - hBot);

            #pragma omp parallel for
                for (int i = 0; i < bottomLayer->getObstacles().size(); ++i){
                    //Check if the layers actually exist as fluid cells, otherwise 0 transport
                    if (topLayer->getObstacles(i) || bottomLayer->getObstacles(i)){
                        middleLayer->setVerticalVelocity(i, vSolid);
                        continue;
                    }

                    //Compute convection intensity

                    //Temperatures
                    V tBot = bottomLayer->getTemperature(i);
                    V tTop = topLayer->getTemperature(i);

                    //Moisture
                    V mBot = bottomLayer->getMoisture(i);
                    V mTop = topLayer->getMoisture(i);

                    //Quantity of air per volume (perfect gas law)
                    V nBot = pBot / universalGasConstant / tBot;
                    V nTop = pTop / universalGasConstant / tTop;

                    //Specific humidities
                    V qBot = mBot / nBot;
                    V qTop = mTop / nTop;

                    //Mean Static Energies
                    V mseBot = Cp * tBot + g * hBot + latVap * qBot;
                    V mseTop = Cp * tTop + g * hTop + latVap * qTop;

                    V vertVel = convectCouplingCoeff * (mseBot - mseTop);
                    vertVel = PWM::Utils::getSign(vertVel) * std::min(debugConvectiveMax, std::abs(vertVel));

                    //Update convectionLayer's velocity.
                    middleLayer->setVerticalVelocity(i, vertVel);
                }
            #pragma omp barrier
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::calcPressureUplift(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer){
            V vSolid = 0.0;//Since mountains don't walk or jump, always 0.
            V avePBot = bottomLayer->getPressure().mean();
            V avePTop = topLayer->getPressure().mean();
            #pragma omp parallel for
                for (int i = 0; i < bottomLayer->getObstacles().size(); ++i){
                    if (topLayer->getObstacles(i) || bottomLayer->getObstacles(i)){
                        middleLayer->setVerticalVelocity(i, vSolid);
                        continue;
                    }
                    V pBot = bottomLayer->getPressure(i) - avePBot;
                    V pTop = topLayer->getPressure(i) - avePTop;
                    V newVV = middleLayer->getVerticalVelocity(i) + (pressureCouplingCoeff * (pTop - pBot));
                    middleLayer->setVerticalVelocity(i, newVV);
                }
            #pragma omp barrier
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::transportFromUplift(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer){
            //Apply vertical transport
            V vSolid = 0.0;//Since mountains don't walk or jump, always 0.
            V heightScaleFactor = topLayer->getHeight() - bottomLayer->getHeight();

            #pragma omp parallel for
                for (int i = 0; i < bottomLayer->getObstacles().size(); ++i){

                    V vertVel = middleLayer->getVerticalVelocity(i);
                    V coeff = (1 - std::exp((-std::abs(vertVel)* getDt()) / heightScaleFactor));

                    V cBot = bottomLayer->getCondensedWater(i);
                    V cTop = topLayer->getCondensedWater(i);
                    V mBot = bottomLayer->getMoisture(i);
                    V mTop = topLayer->getMoisture(i);
                    V pBot = bottomLayer->getPressure(i);
                    V pTop = topLayer->getPressure(i);
                    V tBot = bottomLayer->getTemperature(i);
                    V tTop = topLayer->getTemperature(i);

                    if (topLayer->getObstacles(i) || bottomLayer->getObstacles(i)){
                        bufferPressureB->setData(i, pBot);
                        bufferPressureT->setData(i, pTop);
                        bufferCloudB->setData(i, cBot);
                        bufferCloudT->setData(i, cTop);
                        bufferMoistureB->setData(i, mBot);
                        bufferMoistureT->setData(i, mTop);
                        bufferTemperatureB->setData(i, tBot);
                        bufferTemperatureT->setData(i, tTop);
                        continue;
                    }

                    V quantC = std::abs(cTop - cBot) * coeff;
                    V quantM = std::abs(mTop - mBot) * coeff;
                    V quantT = std::abs(tTop - tBot) * coeff;

                    V newBC, newTC, newBM, newTM, newBT, newTT;
                    if (vertVel > 0.){
                        quantC = std::min(cBot, quantC);
                        newTC = cTop + quantC;
                        newBC = cBot - quantC;
                        bufferCloudB->setData(i, newBC);
                        bufferCloudT->setData(i, newTC);

                        quantM = std::min(mBot, quantM);
                        newTM = mTop + quantM;
                        newBM = mBot - quantM;
                        bufferMoistureB->setData(i, newBM);
                        bufferMoistureT->setData(i, newTM);

                        quantT = std::min(tBot, quantT);
                        newTT = tTop + quantT;
                        newBT = tBot - quantT;
                        bufferTemperatureB->setData(i, newBT);
                        bufferTemperatureT->setData(i, newTT);

                        if (performPressureUplift){
                            V newBP, newTP;
                            newTP = pTop + pressureUpliftCoeff * vertVel;
                            newBP = pBot - pressureUpliftCoeff * vertVel;
                            bufferPressureB->setData(i, newBP);
                            bufferPressureT->setData(i, newTP);
                        }
                        else{
                            bufferPressureB->setData(i, pBot);
                            bufferPressureT->setData(i, pTop);
                        }
                    }
                    else {
                        quantC = std::min(cTop, quantC);
                        newTC = cTop - quantC;
                        newBC = cBot + quantC;
                        bufferCloudB->setData(i, newBC);
                        bufferCloudT->setData(i, newTC);

                        quantM = std::min(mTop, quantM);
                        newTM = mTop - quantM;
                        newBM = mBot + quantM;
                        bufferMoistureB->setData(i, newBM);
                        bufferMoistureT->setData(i, newTM);

                        quantT = std::min(tTop, quantT);
                        newTT = tTop - quantT;
                        newBT = tBot + quantT;
                        bufferTemperatureB->setData(i, newBT);
                        bufferTemperatureT->setData(i, newTT);

                        if (performPressureUplift){
                            V newBP, newTP;
                            newBP = pBot + pressureUpliftCoeff * vertVel;
                            newTP = pTop - pressureUpliftCoeff * vertVel;
                            bufferPressureB->setData(i, newBP);
                            bufferPressureT->setData(i, newTP);
                        }
                        else{
                            bufferPressureB->setData(i, pBot);
                            bufferPressureT->setData(i, pTop);
                        }
                    }
                }
            #pragma omp barrier

            bottomLayer->swapClouds(bufferCloudB);
            bottomLayer->swapMoistures(bufferMoistureB);
            bottomLayer->swapPres(bufferPressureB);
            bottomLayer->swapTemps(bufferTemperatureB);

            topLayer->swapClouds(bufferCloudT);
            topLayer->swapMoistures(bufferMoistureT);
            topLayer->swapPres(bufferPressureT);
            topLayer->swapTemps(bufferTemperatureT);
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::upliftInOutflowApproximation(std::shared_ptr<PWM::Model::airLayer<T, V>>& bottomLayer, std::shared_ptr<PWM::Model::airLayer<T, V>>& topLayer, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& middleLayer){
            std::cerr << "Error! Not implemented for general template, use a specialized template or define upliftInOutflowApproximation(PWM::Model::airLayer<T, V>& bottomLayer, PWM::Model::airLayer<T, V>& topLayer, PWM::Model::convectionLayer<T, V>& middleLayer) for this template!" << std::endl;
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::upliftInOutflowApproximation(){
            std::cerr << "Error! Not implemented for general template, use a specialized template or define upliftInOutflowApproximation() for this template!" << std::endl;
        }

        template<>
        inline void verticalCouplingEngine<PWM::PWMDataStructure::square2DArray<double>, double>::upliftInOutflowApproximation(){
            auto aL = this->getAirLayers();
            auto cL = this->getConvectionLayers();
            int aLCount = aL.size();
            int cellCount = aL[0]->getObstacles().size();
            for (int k = 0; k < aLCount; ++k){
                #pragma omp parallel for
                    for (int i = 0; i < cellCount; ++i){
                        if (aL[k]->getObstacles(i)){
                            continue;
                        }
//                        if (i == 105165){
//                            std::cerr << "Debug point reached" << std::endl;
//                        }
                        double botFlow, topFlow, nettFlow, cPres, newPres;
                        cPres = aL[k]->getPressure(i);
                        if (k == 0){
                            botFlow = 0;
                            topFlow = cL[k]->getVerticalVelocity(i);
                        }else if (k == aLCount - 1){
                            botFlow = cL[k - 1]->getVerticalVelocity(i);
                            topFlow = 0;
                        }else{
                            botFlow = cL[k - 1]->getVerticalVelocity(i);
                            topFlow = cL[k]->getVerticalVelocity(i);
                        }
                        nettFlow = botFlow - topFlow;
                        newPres = cPres + (upliftTransportCoeff * nettFlow);
                        aL[k]->setPressure(i, newPres);
                    }
                #pragma omp barrier
            }
        }

        template<>
        inline void verticalCouplingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::upliftInOutflowApproximation(){
            auto aL = this->getAirLayers();
            auto cL = this->getConvectionLayers();
            int aLCount = aL.size();
            int nX = aL[0]->getObstacles().getX(), nY = aL[0]->getObstacles().getY();
            #pragma omp parallel for
                for (int k = 0; k < aLCount; ++k){
                    for (int i = 0; i < nX; ++i){
                        for (int j = 0; j < nY; ++j){
                            if (aL[k]->getObstacles(i)){
                                continue;
                            }

                            double botFlow, topFlow, nettFlow;
                            if (k == 0){
                                botFlow = 0;
                                topFlow = cL[k]->getVerticalVelocity(i, j);
                            }
                            else if (k == aLCount - 1){
                                botFlow = cL[k - 1]->getVerticalVelocity(i, j);
                                topFlow = 0;
                            }
                            else{
                                botFlow = cL[k - 1]->getVerticalVelocity(i, j);
                                topFlow = cL[k]->getVerticalVelocity(i, j);
                            }
                            nettFlow = botFlow - topFlow;

                            double newVelXTop, newVelXBot, newVelYLeft, newVelYRight;
                            double oldVelXTop, oldVelXBot, oldVelYLeft, oldVelYRight;
                            oldVelXTop = aL[k]->getVelocityTheta(i, j);
                            oldVelXBot = aL[k]->getVelocityTheta(i + 1, j);
                            oldVelYLeft = aL[k]->getVelocityPhi(i, j);
                            oldVelYRight = aL[k]->getVelocityPhi(i, j + 1);
                            newVelXTop = oldVelXTop - (upliftTransportCoeff * nettFlow);
                            newVelXBot = oldVelXBot + (upliftTransportCoeff * nettFlow);
                            newVelYLeft = oldVelYLeft - (upliftTransportCoeff * nettFlow);
                            newVelYRight = oldVelYRight + (upliftTransportCoeff * nettFlow);
                            aL[k]->setVelocityTheta(i, j, newVelXTop);
                            aL[k]->setVelocityTheta(i + 1, j, newVelXBot);
                            aL[k]->setVelocityPhi(i, j, newVelYLeft);
                            aL[k]->setVelocityPhi(i, j + 1, newVelYRight);
                        }
                    }
                }
            #pragma omp barrier
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::smoothVertVelocities(){
            for (auto layer : convectionLayers){
                layer->getVerticalVelocities().movingAverageSmoothing(bufferCloudB, moveAverageWindowRadius);
                layer->swapVertVel(bufferCloudB);
            }
        }

        template<typename T, typename V>
        inline void verticalCouplingEngine<T, V>::step_internal(){
            if (convectionLayers.size() != airLayers.size() - 1){
                std::cerr << "Incorrect layer counts! Must be 1 less convection layer than air layers!" << std::endl;
                return;
            }
            this->sortLayers();
            for (int i = 0; i < convectionLayers.size(); ++i){
                this->calcConvectiveUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
                if (performPressureUplift)
                    this->calcPressureUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
                this->transportFromUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
            }
            this->upliftInOutflowApproximation();
        }

        template<>
        inline void verticalCouplingEngine<PWM::PWMDataStructure::square2DArray<double>, double>::step_internal(){
            if (convectionLayers.size() != airLayers.size() - 1){
                std::cerr << "Incorrect layer counts! Must be 1 less convection layer than air layers!" << std::endl;
                return;
            }
            this->sortLayers();
            for (int i = 0; i < convectionLayers.size(); ++i){
                this->calcConvectiveUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
                if (performPressureUplift)
                    this->calcPressureUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
                this->transportFromUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
            }
            this->upliftInOutflowApproximation();
        }

        template<>
        inline void verticalCouplingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::step_internal(){
            if (convectionLayers.size() != airLayers.size() - 1){
                std::cerr << "Incorrect layer counts! Must be 1 less convection layer than air layers!" << std::endl;
                return;
            }
            this->sortLayers();
            for (int i = 0; i < convectionLayers.size(); ++i){
                this->calcConvectiveUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
                if (performPressureUplift)
                    this->calcPressureUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
                if (performVertVelSmooth)
                    smoothVertVelocities();
                this->transportFromUplift(airLayers[i], airLayers[i + 1], convectionLayers[i]);
            }
            this->upliftInOutflowApproximation();
        }
    }
}

#endif //PWM_ENGINE_VERTICAL_COUPLING_ENGINE_H
