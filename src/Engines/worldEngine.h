#ifndef PWM_WORLD_ENGINE_H
#define PWM_WORLD_ENGINE_H

#include "abstractEngine.h"
#include "advectionEngine.h"
#include "atmosphericHeatingEngine.h"
#include "heatDissipationEngine.h"
#include "moistureUptakeEngine.h"
#include "obstaclesEngine.h"
#include "phaseTransitionEngine.h"
#include "precipitationEngine.h"
#include "pressureEngine.h"
#include "sunEngine.h"
#include "terrainHeatingEngine.h"
#include "verticalCouplingEngine.h"
#include "vizUtils.h"
#include "world.h"

namespace PWM{
    namespace Engine{
        template<typename T, typename TT, typename V, typename VV> class worldEngine : public AbstractEngine{
            private:
                std::shared_ptr<PWM::Model::world<T, TT, V, VV>> worldModel;

                std::shared_ptr<advectionEngine<T, V>> advectEngine;

                std::shared_ptr<atmosphericHeatingEngine<T, TT, V, VV>> atmoHeatEngine;

                std::shared_ptr<heatDissipationEngine<T, V>> heatDissipEngine;

                std::shared_ptr<moistureUptakeEngine<T, TT, V, VV>> moistUptakeEngine;

                std::shared_ptr<obstaclesEngine<T, TT, V, VV>> obsEngine;

                std::shared_ptr<phaseTransitionEngine<T, V>> phaseTransEngine;

                std::shared_ptr<precipitationEngine<T, TT, V, VV>> precipEngine;

                std::shared_ptr<pressureEngine<T, V>> presEngine;

                std::shared_ptr<sunEngine<V>> solarEngine;

                std::shared_ptr<terrainHeatingEngine<T, TT, V, VV>> terrHeatEngine;

                std::shared_ptr<verticalCouplingEngine<T, V>> vertCouplingEngine;

                int stepCount = 0;

                float vertCoupEngDt = 1.0;


                bool checkMoisture() const;
            protected:
                void step_internal() override;
            public:
                worldEngine(float dt = 60, bool active = true);
                worldEngine(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldM, float dt, bool active);
                worldEngine(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldM, float dt, float cveDT, bool active);
                worldEngine(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldM, PWM::Utils::settings& s, bool active);
                int getStepCount() const;

                bool printViz = false;
                std::string outputDirectory = "./";

                const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& getWorldModel() const;

                void printTiming(bool detail) const;
        };

        template<typename T, typename TT, typename V, typename VV>
        inline worldEngine<T, TT, V, VV>::worldEngine(float dt, bool active): AbstractEngine(dt, active), worldModel(nullptr), stepCount(0){
            advectEngine = nullptr;
            atmoHeatEngine = nullptr;
            vertCouplingEngine = nullptr;
            heatDissipEngine = nullptr;
            moistUptakeEngine = nullptr;
            obsEngine = nullptr;
            phaseTransEngine = nullptr;
            precipEngine = nullptr;
            presEngine = nullptr;
            solarEngine = nullptr;
            terrHeatEngine = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline worldEngine<T, TT, V, VV>::worldEngine(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldM, float dt, bool active): AbstractEngine(dt, active), worldModel(worldM), vertCoupEngDt(dt){
            stepCount = 0;
            auto terr = worldModel->getTerrain();
            
            advectEngine = std::make_shared<advectionEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), worldModel->getWorldSize(), dt, active);
            for (auto l : worldModel->getAirLayers())
                advectEngine->addLayer(l);
            
            atmoHeatEngine = std::make_shared<atmosphericHeatingEngine<T, TT, V, VV>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                atmoHeatEngine->addAirLayer(l);
            atmoHeatEngine->setTerrain(terr);
            
            heatDissipEngine = std::make_shared<heatDissipationEngine<T, V>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                heatDissipEngine->addAirLayer(l);
            
            moistUptakeEngine = std::make_shared<moistureUptakeEngine<T, TT, V, VV>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                moistUptakeEngine->addAirLayer(l);
            moistUptakeEngine->setTerrain(terr);
            
            obsEngine = std::make_shared<obstaclesEngine<T, TT, V, VV>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                obsEngine->addAirLayer(l);
            obsEngine->setTerrain(terr);
            
            phaseTransEngine = std::make_shared<phaseTransitionEngine<T, V>>(dt, active, true, true, true, true, 0.001, 0.001);
            for (auto l : worldModel->getAirLayers())
                phaseTransEngine->addAirLayer(l);

            precipEngine = std::make_shared<precipitationEngine<T, TT, V, VV>>(dt, active, 0.05, 0.001, 1.0, 1.0);
            for (auto l : worldModel->getAirLayers())
                precipEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                precipEngine->addConvectionLayer(c);
            precipEngine->setTerrain(terr);

            presEngine = std::make_shared<pressureEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), worldModel->getWorldSize(), dt, active);
            for (auto l : worldModel->getAirLayers())
                presEngine->addLayer(l);


            solarEngine = std::make_shared<sunEngine<V>>(dt, active);
            for (auto s : worldModel->getSuns())
                solarEngine->addSun(s);
            
            terrHeatEngine = std::make_shared<terrainHeatingEngine<T, TT, V, VV>>(dt, active);
            for (auto s : worldModel->getSuns())
                terrHeatEngine->addSun(s);
            terrHeatEngine->setTerrain(terr);

            vertCouplingEngine = std::make_shared<verticalCouplingEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), dt, active);
            for (auto l : worldModel->getAirLayers())
                vertCouplingEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                vertCouplingEngine->addConvectionLayer(c);
        }

        template<>
        inline worldEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::worldEngine(std::shared_ptr<PWM::Model::world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>& worldM, float dt, bool active): AbstractEngine(dt, active), worldModel(worldM), vertCoupEngDt(dt){
            stepCount = 0;
            auto terr = worldModel->getTerrain();
            size_t width = worldModel->getAirLayer(0)->getObstacles().getX();
            size_t height = worldModel->getAirLayer(0)->getObstacles().getY();
            double xSize = worldModel->getAirLayer(0)->getObstacles().getXWorldLength();
            double ySize = worldModel->getAirLayer(0)->getObstacles().getYWorldLength();

            advectEngine = std::make_shared<advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(width, height, xSize, ySize, dt, active, false);
            for (auto l : worldModel->getAirLayers())
                advectEngine->addLayer(l);

            atmoHeatEngine = std::make_shared<atmosphericHeatingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                atmoHeatEngine->addAirLayer(l);
            atmoHeatEngine->setTerrain(terr);

            heatDissipEngine = std::make_shared<heatDissipationEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                heatDissipEngine->addAirLayer(l);

            moistUptakeEngine = std::make_shared<moistureUptakeEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                moistUptakeEngine->addAirLayer(l);
            moistUptakeEngine->setTerrain(terr);

            obsEngine = std::make_shared<obstaclesEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                obsEngine->addAirLayer(l);
            obsEngine->setTerrain(terr);

            phaseTransEngine = std::make_shared<phaseTransitionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(dt, active, true, true, true, true, 0.1, 0.1);
            for (auto l : worldModel->getAirLayers())
                phaseTransEngine->addAirLayer(l);

            precipEngine = std::make_shared<precipitationEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active, 0.05, 0.001, 1.0, 1.0);
            for (auto l : worldModel->getAirLayers())
                precipEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                precipEngine->addConvectionLayer(c);
            precipEngine->setTerrain(terr);

            presEngine = std::make_shared<pressureEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(width, height, xSize, ySize, dt, active);
            for (auto l : worldModel->getAirLayers())
                presEngine->addLayer(l);

            solarEngine = std::make_shared<sunEngine<double>>(dt, active);
            for (auto s : worldModel->getSuns())
                solarEngine->addSun(s);

            terrHeatEngine = std::make_shared<terrainHeatingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto s : worldModel->getSuns())
                terrHeatEngine->addSun(s);
            terrHeatEngine->setTerrain(terr);

            vertCouplingEngine = std::make_shared<verticalCouplingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(width, height, xSize, ySize, dt, active, false);
            // JG last flag switches off pressure uplift
            for (auto l : worldModel->getAirLayers())
                vertCouplingEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                vertCouplingEngine->addConvectionLayer(c);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline worldEngine<T, TT, V, VV>::worldEngine(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldM, float dt, float cveDT, bool active): AbstractEngine(dt, active), worldModel(worldM), vertCoupEngDt(cveDT){
            stepCount = 0;
            auto terr = worldModel->getTerrain();

            advectEngine = std::make_shared<advectionEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), worldModel->getWorldSize(), dt, active);
            for (auto l : worldModel->getAirLayers())
                advectEngine->addLayer(l);

            atmoHeatEngine = std::make_shared<atmosphericHeatingEngine<T, TT, V, VV>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                atmoHeatEngine->addAirLayer(l);
            atmoHeatEngine->setTerrain(terr);

            heatDissipEngine = std::make_shared<heatDissipationEngine<T, V>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                heatDissipEngine->addAirLayer(l);

            moistUptakeEngine = std::make_shared<moistureUptakeEngine<T, TT, V, VV>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                moistUptakeEngine->addAirLayer(l);
            moistUptakeEngine->setTerrain(terr);

            obsEngine = std::make_shared<obstaclesEngine<T, TT, V, VV>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                obsEngine->addAirLayer(l);
            obsEngine->setTerrain(terr);

            phaseTransEngine = std::make_shared<phaseTransitionEngine<T, V>>(dt, active, true, true, true, true, 0.001, 0.001);
            for (auto l : worldModel->getAirLayers())
                phaseTransEngine->addAirLayer(l);

            precipEngine = std::make_shared<precipitationEngine<T, TT, V, VV>>(dt, active, 0.05, 0.001, 1.0, 1.0);
            for (auto l : worldModel->getAirLayers())
                precipEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                precipEngine->addConvectionLayer(c);
            precipEngine->setTerrain(terr);

            presEngine = std::make_shared<pressureEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), worldModel->getWorldSize(), dt, active);
            for (auto l : worldModel->getAirLayers())
                presEngine->addLayer(l);

            solarEngine = std::make_shared<sunEngine<V>>(dt, active);
            for (auto s : worldModel->getSuns())
                solarEngine->addSun(s);

            terrHeatEngine = std::make_shared<terrainHeatingEngine<T, TT, V, VV>>(dt, active);
            for (auto s : worldModel->getSuns())
                terrHeatEngine->addSun(s);
            terrHeatEngine->setTerrain(terr);

            vertCouplingEngine = std::make_shared<verticalCouplingEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), vertCoupEngDt, active);
            for (auto l : worldModel->getAirLayers())
                vertCouplingEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                vertCouplingEngine->addConvectionLayer(c);
        }

        template<>
        inline worldEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::worldEngine(std::shared_ptr<PWM::Model::world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>& worldM, float dt, float cveDT, bool active): AbstractEngine(dt, active), worldModel(worldM), vertCoupEngDt(cveDT){
            stepCount = 0;
            auto terr = worldModel->getTerrain();
            size_t width = worldModel->getAirLayer(0)->getObstacles().getX();
            size_t height = worldModel->getAirLayer(0)->getObstacles().getY();
            double xSize = worldModel->getAirLayer(0)->getObstacles().getXWorldLength();
            double ySize = worldModel->getAirLayer(0)->getObstacles().getYWorldLength();

            advectEngine = std::make_shared<advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(width, height, xSize, ySize, dt, active, false);
            for (auto l : worldModel->getAirLayers())
                advectEngine->addLayer(l);

            atmoHeatEngine = std::make_shared<atmosphericHeatingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                atmoHeatEngine->addAirLayer(l);
            atmoHeatEngine->setTerrain(terr);

            heatDissipEngine = std::make_shared<heatDissipationEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                heatDissipEngine->addAirLayer(l);

            moistUptakeEngine = std::make_shared<moistureUptakeEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                moistUptakeEngine->addAirLayer(l);
            moistUptakeEngine->setTerrain(terr);

            obsEngine = std::make_shared<obstaclesEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto l : worldModel->getAirLayers())
                obsEngine->addAirLayer(l);
            obsEngine->setTerrain(terr);

            phaseTransEngine = std::make_shared<phaseTransitionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(dt, active, true, true, true, true, 0.001, 0.001);
            for (auto l : worldModel->getAirLayers())
                phaseTransEngine->addAirLayer(l);

            precipEngine = std::make_shared<precipitationEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active, 0.05, 0.001, 1.0, 1.0);
            for (auto l : worldModel->getAirLayers())
                precipEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                precipEngine->addConvectionLayer(c);
            precipEngine->setTerrain(terr);

            presEngine = std::make_shared<pressureEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(width, height, xSize, ySize, dt, active);
            for (auto l : worldModel->getAirLayers())
                presEngine->addLayer(l);

            solarEngine = std::make_shared<sunEngine<double>>(dt, active);
            for (auto s : worldModel->getSuns())
                solarEngine->addSun(s);

            terrHeatEngine = std::make_shared<terrainHeatingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(dt, active);
            for (auto s : worldModel->getSuns())
                terrHeatEngine->addSun(s);
            terrHeatEngine->setTerrain(terr);

            vertCouplingEngine = std::make_shared<verticalCouplingEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(width, height, xSize, ySize, vertCoupEngDt, active, false);
            for (auto l : worldModel->getAirLayers())
                vertCouplingEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                vertCouplingEngine->addConvectionLayer(c);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline worldEngine<T, TT, V, VV>::worldEngine(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldM, PWM::Utils::settings& s, bool active): AbstractEngine(s.weDt, active), worldModel(worldM){
            stepCount = 0;
            auto terr = worldModel->getTerrain();

            advectEngine = std::make_shared<advectionEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), worldModel->getWorldSize(), s.weDt, active);
            for (auto l : worldModel->getAirLayers())
                advectEngine->addLayer(l);
            advectEngine->setCoriolis(s.coriolis);

            atmoHeatEngine = std::make_shared<atmosphericHeatingEngine<T, TT, V, VV>>(s.weDt, active);
            for (auto l : worldModel->getAirLayers())
                atmoHeatEngine->addAirLayer(l);
            atmoHeatEngine->setTerrain(terr);
            atmoHeatEngine->setSurfaceFluxRatio(s.aheSurfaceFluxRatio);
            atmoHeatEngine->setThreshold(s.aheThreshold);

            vertCouplingEngine = std::make_shared<verticalCouplingEngine<T, V>>(worldModel->getAirLayer(0)->getObstacles().getWidth(), s.weDt, active);
            for (auto l : worldModel->getAirLayers())
                vertCouplingEngine->addAirLayer(l);
            for (auto c : worldModel->getConvectionLayers())
                vertCouplingEngine->addConvectionLayer(c);
            vertCouplingEngine->setConvectCouplingCoeff(s.vceConvectCouplingCoeff);
            vertCouplingEngine->setPressureCouplingCoeff(s.vcePressureCouplingCoeff);
            vertCouplingEngine->setPressureUpliftCoeff(s.vcePressureUpliftCoeff);
            vertCouplingEngine->setPressureUplift(s.vcePerformPressureUplift);

            heatDissipEngine = std::make_shared<heatDissipationEngine<T, V>>(s.weDt, active);
            for (auto l : worldModel->getAirLayers())
                heatDissipEngine->addAirLayer(l);
            heatDissipEngine->setCoefficient(s.hdeCoefficient);

            moistUptakeEngine = std::make_shared<moistureUptakeEngine<T, TT, V, VV>>(s.weDt, active);
            for (auto l : worldModel->getAirLayers())
                moistUptakeEngine->addAirLayer(l);
            moistUptakeEngine->setTerrain(terr);

            obsEngine = std::make_shared<obstaclesEngine<T, TT, V, VV>>(s.weDt, active);
            for (auto l : worldModel->getAirLayers())
                obsEngine->addAirLayer(l);
            obsEngine->setTerrain(terr);

            phaseTransEngine = std::make_shared<phaseTransitionEngine<T, V>>(s.weDt, active,
                                                                             s.pteDoEvaporation,
                                                                             s.pteDoCondensation,
                                                                             s.pteDoEvaporationHeatExchange,
                                                                             s.pteDoCondensationHeatExchange,
                                                                             s.pteEvaporationCoefficient,
                                                                             s.pteCondensationCoefficient);
            for (auto l : worldModel->getAirLayers())
                phaseTransEngine->addAirLayer(l);

            solarEngine = std::make_shared<sunEngine<V>>(s.weDt, active);
            for (auto s : worldModel->getSuns())
                solarEngine->addSun(s);

            terrHeatEngine = std::make_shared<terrainHeatingEngine<T, TT, V, VV>>(s.weDt, active);
            for (auto s : worldModel->getSuns())
                terrHeatEngine->addSun(s);
            terrHeatEngine->setTerrain(terr);
            terrHeatEngine->setHeatLossCoefficient(s.theHeatLossCoefficient);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline bool worldEngine<T, TT, V, VV>::checkMoisture() const{
            auto layers = this->getWorldModel()->getAirLayers();
            for (auto i : layers){
                V min = i->getMoisture().min();
                if (min < 0)
                    return false;
            }
            return true;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void worldEngine<T, TT, V, VV>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define worldEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void worldEngine<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::step_internal(){
            //Sun moves
            solarEngine->step();

            //Sun heats terrain
            terrHeatEngine->step();
            if (printViz){
                std::stringstream f;
                    f << outputDirectory << "/Terrain_Temp_During_Step_" << stepCount << "_After_TerrHeatEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f.str(), worldModel->getTerrain()->getTemperature());
            }

            //Ground heats the lowest air layer
            atmoHeatEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f;
                    f << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_AtmoHeatEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f.str(), worldModel->getAirLayer(i)->getTemperature());
                }
            }

            //Ground moisture evaporates into air
            moistUptakeEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f;
                    f << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_MoistUptakeEngineStep.ppm";
                    PWM::Utils::writeMoisImage<double>(f.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getMoisture());
                }
            }

            //Air moisture and particulates drop out of the air
            precipEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_Ash_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Water_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                    PWM::Utils::writeMoisImage<double>(f1.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getParticulates());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getCondensedWater());
                    if (i < worldModel->getAirLayers().size() - 1){
                        std::stringstream f3, f4, f5;
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_VertVelocity_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_Ashfall_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_Rainfall_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getVerticalVelocities());
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getRainfall());
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getAshfall());
                    }
                }
            }

            //Heat dissipates out (because entropy) - an approximation of radiative cooling
//            heatDissipEngine->step();

            //Magic of fluid dynamics! Advection happens in each layer
            advectEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2, f3, f4;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f3 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f4 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getMoisture());
                    PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage<double>(f4.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            //Obstacles kill waves, so do that.
            obsEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_ObsEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_ObsEngineStep.ppm";
                    PWM::Utils::writeVelImage(f1.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage(f2.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            //Fix the divergence that got broken
            presEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2, f3, f4;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    f3 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    f4 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getMoisture());
                    PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage<double>(f4.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            //Layers have to talk to each other, and convection is the line.
            vertCouplingEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_ConvectEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_ConvectEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getMoisture());
                    if (i < worldModel->getAirLayers().size() - 1){
                        std::stringstream f3;
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_VertVelocity_During_Step_" << stepCount << "_After_ConvectEngineStep.ppm";
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getVerticalVelocities());
                    }
                }
            }

            //Fix values broken because Navier-Stokes is a *****...
            presEngine->smoothPressure();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_PressureSmoothingStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_PressureSmoothingStep.ppm";
                    PWM::Utils::writeVelImage(f1.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage(f2.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            //Now apply the physics to the air just moved around, phase transitions
            phaseTransEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getPlanet(), worldModel->getAirLayer(i)->getMoisture());
                }
            }

            ++stepCount;
        }

        template<>
        inline void worldEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::step_internal(){
            //Sun moves
            solarEngine->step();

            //Sun heats terrain
            terrHeatEngine->step();
            if (printViz){
                std::stringstream f;
                    f << outputDirectory << "/Terrain_Temp_During_Step_" << stepCount << "_After_TerrHeatEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f.str(), worldModel->getTerrain()->getTemperature());
            }

            //Ground heats the lowest air layer
            atmoHeatEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f;
                    f << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_AtmoHeatEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f.str(), worldModel->getAirLayer(i)->getTemperature());
                }
            }

            //Ground moisture evaporates into air
            moistUptakeEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f;
                    f << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_MoistUptakeEngineStep.ppm";
                    PWM::Utils::writeMoisImage<double>(f.str(), worldModel->getAirLayer(i)->getMoisture());
                }
            }

            //Air moisture and particulates drop out of the air
            precipEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_Ash_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Water_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                    PWM::Utils::writeMoisImage<double>(f1.str(), worldModel->getAirLayer(i)->getParticulates());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getAirLayer(i)->getCondensedWater());
                    if (i < worldModel->getAirLayers().size() - 1){
                        std::stringstream f3, f4, f5;
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_VertVelocity_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_Ashfall_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                        f3 << outputDirectory << "/ConvectionLayer_" << i << "_Rainfall_During_Step_" << stepCount << "_After_PrecipEngineStep.ppm";
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getVerticalVelocities());
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getRainfall());
                        PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getAshfall());
                    }
                }
            }

            //Heat dissipates out (because entropy) - an approximation of radiative cooling
//            heatDissipEngine->step();

            //Magic of fluid dynamics! Advection happens in each layer
            advectEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2, f3, f4;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f3 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    f4 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_AdvectEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getAirLayer(i)->getMoisture());
                    PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage<double>(f4.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            //Obstacles kill waves, so do that.
            obsEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2;
                    f1 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_ObsEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_ObsEngineStep.ppm";
                    PWM::Utils::writeVelImage(f1.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage(f2.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            //Fix the divergence that got broken
            presEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2, f3, f4;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    f3 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    f4 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_DivergeCorrectStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getAirLayer(i)->getMoisture());
                    PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                    PWM::Utils::writeVelImage<double>(f4.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                }
            }

            if (! std::fmod(getSimTimePassed(), vertCoupEngDt)){
                //Layers have to talk to each other, and convection is the line.
                vertCouplingEngine->step();
                if (printViz){
                    for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                        std::stringstream f1, f2;
                        f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_ConvectEngineStep.ppm";
                        f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_ConvectEngineStep.ppm";
                        PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                        PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getAirLayer(i)->getMoisture());
                        if (i < worldModel->getAirLayers().size() - 1){
                            std::stringstream f3;
                            f3 << outputDirectory << "/ConvectionLayer_" << i << "_VertVelocity_During_Step_" << stepCount << "_After_ConvectEngineStep.ppm";
                            PWM::Utils::writeVelImage<double>(f3.str(), worldModel->getConvectionLayer(i)->getVerticalVelocities());
                        }
                    }
                }

                //Fix values broken because Navier-Stokes is a *****...
//                presEngine->solvePressureProjection();
                if (printViz){
                    for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                        std::stringstream f1, f2;
                        f1 << outputDirectory << "/Layer_" << i << "_VelocityX_During_Step_" << stepCount << "_After_PressureSmoothingStep.ppm";
                        f2 << outputDirectory << "/Layer_" << i << "_VelocityY_During_Step_" << stepCount << "_After_PressureSmoothingStep.ppm";
                        PWM::Utils::writeVelImage(f1.str(), worldModel->getAirLayer(i)->getVelocityPhi());
                        PWM::Utils::writeVelImage(f2.str(), worldModel->getAirLayer(i)->getVelocityTheta());
                    }
                }
            }

            //Now apply the physics to the air just moved around, phase transitions
            phaseTransEngine->step();
            if (printViz){
                for (int i = 0; i < worldModel->getAirLayers().size(); ++i){
                    std::stringstream f1, f2, f3;
                    f1 << outputDirectory << "/Layer_" << i << "_Temp_During_Step_" << stepCount << "_After_PhaseTransitionEngineStep.ppm";
                    f2 << outputDirectory << "/Layer_" << i << "_Moisture_During_Step_" << stepCount << "_After_PhaseTransitionEngineStep.ppm";
                    f3 << outputDirectory << "/Layer_" << i << "_Condensed_Water_During_Step_" << stepCount << "_After_PhaseTransitionEngineStep.ppm";
                    PWM::Utils::writeTempImage<double>(f1.str(), worldModel->getAirLayer(i)->getTemperature());
                    PWM::Utils::writeMoisImage<double>(f2.str(), worldModel->getAirLayer(i)->getMoisture());
                    PWM::Utils::writeMoisImage(f3.str(), worldModel->getAirLayer(i)->getCondensedWater());
                }
            }

            ++stepCount;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& worldEngine<T, TT, V, VV>::getWorldModel() const{
            return worldModel;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline int worldEngine<T, TT, V, VV>::getStepCount() const{
            return stepCount;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void worldEngine<T, TT, V, VV>::printTiming(bool detail) const{
            std::cout << "Error! Not implemented for general template, use a specialized template or define worldEngine::printTiming(bool detail) for this template!" << std::endl;
        }

        template<>
        inline void worldEngine<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::printTiming(bool detail) const{
            auto width = getWorldModel()->getWorldCellWidth();
            double worldTime = getRunTimePassed();
            std::cout << "Time taken for " << getStepCount() << " steps of WorldEngine with a simulation timestep of " << getDt() << " seconds and " << worldModel->getAirLayers().size() << " layers on a " << width << "x" << width << " grid: " << worldTime << " seconds." << std::endl;
            if (detail){
                double advectTime = advectEngine->getRunTimePassed();
                std::cout << "Advection time:\t\t\t" << advectTime << " s\t" << (advectTime / worldTime) * 100 << "%" << std::endl;
                double atmoHeatTime = atmoHeatEngine->getRunTimePassed();
                std::cout << "Atmospheric heating time:\t" << atmoHeatTime << " s\t" << (atmoHeatTime / worldTime) * 100 << "%" << std::endl;
                double heatDissipTime = heatDissipEngine->getRunTimePassed();
                std::cout << "Heat dissipation time:\t\t" << heatDissipTime << " s\t" << (heatDissipTime / worldTime) * 100 << "%" << std::endl;
                double moistUptakeTime = moistUptakeEngine->getRunTimePassed();
                std::cout << "Moisture uptake time:\t\t" << moistUptakeTime << " s\t" << (moistUptakeTime / worldTime) * 100 << "%" << std::endl;
                double obsTime = obsEngine->getRunTimePassed();
                std::cout << "Obstacles time:\t\t\t" << obsTime << " s\t" << (obsTime / worldTime) * 100 << "%" << std::endl;
                double phaseTransTime = phaseTransEngine->getRunTimePassed();
                std::cout << "Phase transition time:\t\t" << phaseTransTime << " s\t" << (phaseTransTime / worldTime) * 100 << "%" << std::endl;
                double precipTime = precipEngine->getRunTimePassed();
                std::cout << "Phase transition time:\t\t" << precipTime << " s\t" << (precipTime / worldTime) * 100 << "%" << std::endl;
                double presTime = presEngine->getRunTimePassed();
                std::cout << "Pressure solve time:\t\t" << presTime << " s\t" << (presTime / worldTime) * 100 << "%" << std::endl;
                double solarTime = solarEngine->getRunTimePassed();
                std::cout << "Solar cycle time:\t\t" << solarTime << " s\t" << (solarTime / worldTime) * 100 << "%" << std::endl;
                double terrHeatTime = terrHeatEngine->getRunTimePassed();
                std::cout << "Terrain heating time:\t\t" << terrHeatTime << " s\t" << (terrHeatTime / worldTime) * 100 << "%" << std::endl;
                double vertCouplTime = vertCouplingEngine->getRunTimePassed();
                std::cout << "Vertical coupling time:\t\t" << vertCouplTime << " s\t" << (vertCouplTime / worldTime) * 100 << "%" << std::endl;
            }
        }

        template<>
        inline void worldEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::printTiming(bool detail) const{
            auto widthX = getWorldModel()->getAirLayer(0)->getObstacles().getX();
            auto widthY = getWorldModel()->getAirLayer(0)->getObstacles().getY();
            double worldTime = getRunTimePassed();
            std::cout << "Time taken for " << getStepCount() << " steps of WorldEngine with a simulation timestep of " << getDt() << " seconds and " << worldModel->getAirLayers().size() << " layers on a " << widthX << "x" << widthY << " grid: " << worldTime << " seconds." << std::endl;
            if (detail){
                double advectTime = advectEngine->getRunTimePassed();
                std::cout << "Advection time:\t\t\t" << advectTime << " s\t" << (advectTime / worldTime) * 100 << "%" << std::endl;
                double atmoHeatTime = atmoHeatEngine->getRunTimePassed();
                std::cout << "Atmospheric heating time:\t" << atmoHeatTime << " s\t" << (atmoHeatTime / worldTime) * 100 << "%" << std::endl;
                double heatDissipTime = heatDissipEngine->getRunTimePassed();
                std::cout << "Heat dissipation time:\t\t" << heatDissipTime << " s\t" << (heatDissipTime / worldTime) * 100 << "%" << std::endl;
                double moistUptakeTime = moistUptakeEngine->getRunTimePassed();
                std::cout << "Moisture uptake time:\t\t" << moistUptakeTime << " s\t" << (moistUptakeTime / worldTime) * 100 << "%" << std::endl;
                double obsTime = obsEngine->getRunTimePassed();
                std::cout << "Obstacles time:\t\t\t" << obsTime << " s\t" << (obsTime / worldTime) * 100 << "%" << std::endl;
                double phaseTransTime = phaseTransEngine->getRunTimePassed();
                std::cout << "Phase transition time:\t\t" << phaseTransTime << " s\t" << (phaseTransTime / worldTime) * 100 << "%" << std::endl;
                double precipTime = precipEngine->getRunTimePassed();
                std::cout << "Precipitation time:\t\t" << precipTime << " s\t" << (precipTime / worldTime) * 100 << "%" << std::endl;
                double presTime = presEngine->getRunTimePassed();
                std::cout << "Pressure solve time:\t\t" << presTime << " s\t" << (presTime / worldTime) * 100 << "%" << std::endl;
                double solarTime = solarEngine->getRunTimePassed();
                std::cout << "Solar cycle time:\t\t" << solarTime << " s\t" << (solarTime / worldTime) * 100 << "%" << std::endl;
                double terrHeatTime = terrHeatEngine->getRunTimePassed();
                std::cout << "Terrain heating time:\t\t" << terrHeatTime << " s\t" << (terrHeatTime / worldTime) * 100 << "%" << std::endl;
                double vertCouplTime = vertCouplingEngine->getRunTimePassed();
                std::cout << "Vertical coupling time:\t\t" << vertCouplTime << " s\t" << (vertCouplTime / worldTime) * 100 << "%" << std::endl;
            }
        }
    }
}

#endif //PWM_WORLD_ENGINE_H
