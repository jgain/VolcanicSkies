#ifndef PWM_PHASE_TRANSITION_ENGINE_H
#define PWM_PHASE_TRANSITION_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include "atmofuncs.h"
#include <omp.h>

namespace PWM{
    namespace Engine{
        template<typename T, typename V> class phaseTransitionEngine : public AbstractEngine{
            private:
                //Toggles on what parts of this engine should be run
                bool doEvaporation;
                bool doCondensation;
                bool doEvaporationHeatExchange;
                bool doCondensationHeatExchange;
                bool doParticleCondensation;


                //Coefficients that affect this engine's calculations.
                V evaporationCoefficient = 0.1f;
                V condensationCoefficient = 0.1f;

                //The vector of air layers that this engine applies to.
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;
            protected:
                void step_internal();
            public:
                phaseTransitionEngine(float dt, bool active);
                phaseTransitionEngine(float dt, bool active, bool doEvap, bool doCondense, bool doEvapHX, bool doCondHX, V evapCoeff, V condCoeff);

                void setDoEvaporation(bool evap);
                void setDoCondensation(bool cond);
                void setDoEvaporationHX(bool evapHX);
                void setDoCondensationHX(bool condHX);
                void setDoParticleCondensation(bool partCond);

                void setEvaporationCoefficient(V evapCoeff);
                void setCondensationCoefficient(V condCoeff);

                bool getDoEvaporation() const;
                bool getDoCondensation() const;
                bool getDoEvaporationHX() const;
                bool getDoCondensationHX() const;
                bool getDoParticleCondensation() const;

                const V getEvaporationCoefficient() const;
                const V getCondensationCoefficient() const;

                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l);
        };

        template<typename T, typename V>
        inline phaseTransitionEngine<T, V>::phaseTransitionEngine(float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
        }

        template<typename T, typename V>
        inline phaseTransitionEngine<T, V>::phaseTransitionEngine(float dt, bool active, bool doEvap, bool doCondense, bool doEvapHX, bool doCondHX, V evapCoeff, V condCoeff) :
                    AbstractEngine(dt, active),
                    doEvaporation(doEvap), doCondensation(doCondense),
                    doEvaporationHeatExchange(doEvapHX), doCondensationHeatExchange(doCondHX),
                    evaporationCoefficient(evapCoeff), condensationCoefficient(condCoeff){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l){
            airLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& phaseTransitionEngine<T, V>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename V>
        inline bool phaseTransitionEngine<T, V>::getDoEvaporation() const{
            return doEvaporation;
        }

        template<typename T, typename V>
        inline bool phaseTransitionEngine<T, V>::getDoCondensation() const{
            return doCondensation;
        }

        template<typename T, typename V>
        inline bool phaseTransitionEngine<T, V>::getDoEvaporationHX() const{
            return doEvaporationHeatExchange;
        }

        template<typename T, typename V>
        inline bool phaseTransitionEngine<T, V>::getDoCondensationHX() const{
            return doCondensationHeatExchange;
        }

        template<typename T, typename V>
        inline bool phaseTransitionEngine<T, V>::getDoParticleCondensation() const{
            return doParticleCondensation;
        }

        template<typename T, typename V>
        inline const V phaseTransitionEngine<T, V>::getEvaporationCoefficient() const{
            return evaporationCoefficient;
        }

        template<typename T, typename V>
        inline const V phaseTransitionEngine<T, V>::getCondensationCoefficient() const{
            return condensationCoefficient;
        }
        
        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setDoEvaporation(bool evap){
            doEvaporation = evap;
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setDoCondensation(bool cond){
            doCondensation = cond;
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setDoEvaporationHX(bool evapHX){
            doEvaporationHeatExchange = evapHX;
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setDoCondensationHX(bool condHX){
            doCondensationHeatExchange = condHX;
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setDoParticleCondensation(bool partCond){
            doParticleCondensation = partCond;
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setEvaporationCoefficient(V evapCoeff){
            evaporationCoefficient = evapCoeff;
        }

        template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::setCondensationCoefficient(V condCoeff){
            condensationCoefficient = condCoeff;
        }

       template<typename T, typename V>
        inline void phaseTransitionEngine<T, V>::step_internal(){
            for (auto l : airLayers){
                V altitude = l->getHeight();
                V atmPres = PWM::Utils::altitudeAdjustedPressure(altitude, l->getPlanet());
                V latHeatVap = l->getPlanet()->getLatentHeatofVaporisation();
                V airHeatCap = l->getPlanet()->getAirHeatCapacity();
                V moistBoilTemp = l->getPlanet()->getMoistureBoilTemp();

                #pragma omp parallel for
                    for (int i = 0; i < l->getObstacles().size(); ++i){
                        if (l->getObstacles(i)){
                            l->setMoisture(i, 0);
                            l->setCondensedWater(i, 0);
                            continue;
                        }
//                        if (i == 128661)
//                            std::cout << "Debug point reached" << std::endl;
                        V temp;
                        V actualTemp;
                        temp = actualTemp = l->getTemperature().getData(i);

                        V mois = l->getMoisture().getData(i);
                        V cond = l->getCondensedWater().getData(i);

                        V satVapPres = PWM::Utils::saturationVapourPressure(temp);

                        V density = PWM::Utils::pressureAdjustedDensity(atmPres, actualTemp, l->getPlanet());
                        V humidity = mois / density;

                        V partialPressure = atmPres * (humidity / (0.622 + humidity));
                        V satHumidity = (0.622 * satVapPres) / (atmPres - satVapPres);
//                        if (temp > 373.15)
//                            std::cout << "Reached debug point at cell " << i << "." << std::endl;

                        //Check if things should be boiling, e.g., middle of volcano plume
                        if ((cond > 0) && (atmPres > partialPressure) && (satVapPres > atmPres)){
                            if (doEvaporation){
                                V water = cond;
                                V newMoisVal = mois + water;
                                l->getMoisture().setData(i, newMoisVal);
                                V newCondVal = 0;
                                l->getCondensedWater().setData(i, newCondVal);
                                l->setClouds(i, newCondVal);

                                if (doEvaporationHeatExchange){
                                    //std::cerr << "Evaporation heat exchange taking place" << std::endl;
                                    V dH = latHeatVap * water;
                                    V deltaT = -dH / (airHeatCap * density);
                                    V newTempVal = actualTemp + deltaT;
                                    l->getTemperature().setData(i, newTempVal);
                                }
                            }
                        }
                        //check if water vapour should condense instead
                        else if ((humidity > 0) && (satVapPres < partialPressure) && (atmPres > satVapPres)){
                            if (doCondensation){
                                V water, diff = humidity - satHumidity;

                                if (doParticleCondensation)
                                    water = diff * condensationCoefficient * (1 + l->getParticulates(i));
                                else
                                    water = diff * condensationCoefficient;

                                if (mois < water)
                                    water = mois;
                                V newMoisVal = mois - water;
                                l->getMoisture().setData(i, newMoisVal);
                                V newCondVal = cond + water;
                                l->getCondensedWater().setData(i, newCondVal);

                                if (doCondensationHeatExchange){
                                    //std::cerr << "Condensation heat exchange taking place" << std::endl;
                                    V dH = latHeatVap * water;
                                    V deltaT = dH / (airHeatCap * density);
                                    V newTempVal = actualTemp + deltaT;
                                    l->getTemperature().setData(i, newTempVal);
                                }
                            }
                        }
                        //check if water should evaporate
                        else if ((cond > 0) && (satVapPres > partialPressure) && (atmPres > satVapPres)){
                            if (doEvaporation){
                                V diff = satHumidity - humidity;
                                V water = diff * evaporationCoefficient;
                                if (cond < water)
                                    water = cond;
                                V newMoisVal = mois + water;
                                l->getMoisture().setData(i, newMoisVal);
                                V newCondVal = cond - water;
                                l->getCondensedWater().setData(i, newCondVal);

                                if (doEvaporationHeatExchange){
                                    //std::cerr << "Evaporation heat exchange taking place" << std::endl;
                                    V dH = latHeatVap * water;
                                    V deltaT = -dH / (airHeatCap * density);
                                    V newTempVal = actualTemp + deltaT;
                                    l->getTemperature().setData(i, newTempVal);
                                }
                            }
                        }
                    }
                #pragma omp barrier
            }
        }
    }
}

#endif //PWM_PHASE_TRANSITION_ENGINE_H
