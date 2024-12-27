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
#ifndef PWM_MODEL_TERRAIN_H
#define PWM_MODEL_TERRAIN_H

#include "flatStaggeredGrid.h"
#include "planet.h"
#include "square2DArray.h"

namespace PWM{
    namespace Model{
        template<typename T, typename TT, typename V, typename VV> class terrain{
            private:
                //The planet model for this terrain
                std::shared_ptr<planet> Planet;

                //The heightmap of the terrain, used for rendering and
                //allocating the 'Obstacles' part of the airlayers. In metres.
                std::shared_ptr<T> elevation;

                //The temperature of the terrain in Kelvin
                std::shared_ptr<T> temperature;

                //The moisture of the terrain in kg per cubic metre
                std::shared_ptr<T> moisture;

                //The ash on top of the terrain in kg per cubic metre
                std::shared_ptr<T> ash;

                //The terrain type at each cell. This determines
                //the albedo and other factors that influence
                //weather.
                std::shared_ptr<TT> cellTerrainType;
            
            public:
                //Constructors
                //Default constructor not meant to be used
                terrain(std::shared_ptr<planet> p);

                //Constructor meant to be used with square2DArray.h
                terrain(std::shared_ptr<planet> p, const size_t width, const V size = 50000.f);


                //Constructor meant to be used with square2DArray.h
                terrain(std::shared_ptr<planet> p, const size_t xWidth, const size_t yHeight, const V xSize, const V ySize);

                //Copy constructor and copy assignment
                terrain(terrain<T, TT, V, VV>& other);
                terrain<T, TT, V, VV>& operator=(terrain<T, TT, V, VV>& other);

                const T& getElevation() const;
                const T& getTemperature() const;
                const T& getMoisture() const;
                const T& getAsh() const;
                const TT& getTerrainType() const;
                const std::shared_ptr<planet> getPlanet() const;

                const V getElevation(int index) const;
                const V getTemperature(int index) const;
                const V getMoisture(int index) const;
                const V getAsh(int index) const;
                const VV getTerrainType(int index) const;

                //For modifications from the phasechange and vegetation engines
                void setTemperature(int index, V& newTemp);
                void setMoisture(int index, V& newMoist);
                void setAsh(int index, V& newAsh);
                void addAsh(int index, V& ashAdded);
                void setTerrainType(int index, VV& newType);
                void setElevation(int index, V& newElev);
                
                void init(std::vector<std::pair<V, bool>>& terrainData);

                bool operator==(terrain<T, TT, V, VV>& other) const;
                bool operator!=(terrain<T, TT, V, VV>& other) const;
        };

        template<typename T, typename TT, typename V, typename VV>
        inline terrain<T, TT, V, VV>::terrain(std::shared_ptr<planet> p){
            elevation = nullptr;
            temperature = nullptr;
            moisture = nullptr;
            ash = nullptr;
            cellTerrainType = nullptr;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline terrain<T, TT, V, VV>::terrain(std::shared_ptr<planet> p, const size_t width, const V size) : Planet(p){
            elevation = std::make_shared<T>(T(width, size));
            temperature = std::make_shared<T>(T(width, size));
            moisture = std::make_shared<T>(T(width, size));
            ash = std::make_shared<T>(T(width, size));
            cellTerrainType = std::make_shared<TT>(TT(width, size));
            temperature->copy(Planet->getAveTemp());
            V x = V();
            elevation->copy(x);
            moisture->copy(x);
            ash->copy(x);
            VV y = VV();
            cellTerrainType->copy(y);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline terrain<T, TT, V, VV>::terrain(std::shared_ptr<planet> p, const size_t xWidth, const size_t yHeight, const V xSize, const V ySize) : Planet(p){
            elevation = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            temperature = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            moisture = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            ash = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            cellTerrainType = std::make_shared<TT>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            temperature->copy(Planet->getAveTemp());
            V x = V();
            elevation->copy(x);
            moisture->copy(x);
            VV y = VV();
            cellTerrainType->copy(y);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline terrain<T, TT, V, VV>::terrain(terrain<T, TT, V, VV>& other) : Planet(other.getPlanet()){
            size_t width = other.getElevation().getWidth();
            V size = other.getElevation().getSize();
            elevation = std::make_shared<T>(T(width, size));
            temperature = std::make_shared<T>(T(width, size));
            moisture = std::make_shared<T>(T(width, size));
            ash = std::make_shared<T>(T(width, size));
            cellTerrainType = std::make_shared<TT>(TT(width, size));
            
            elevation->copy(other.getElevation());
            temperature->copy(other.getTemperature());
            moisture->copy(other.getMoisture());
            ash->copy(other.getAsh());
            cellTerrainType->copy(other.getTerrainType());
        }

        template<>
        inline terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::terrain(terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>& other){
            this->Planet = other.getPlanet();
            size_t nX = other.getElevation().getX();
            size_t nY = other.getElevation().getY();
            double xSize = other.getElevation().getXWorldLength();
            double ySize = other.getElevation().getYWorldLength();
            elevation = std::make_shared<PWM::PWMDataStructure::flatStaggeredGrid<double>>(nX, nY, 0.5, 0.5, xSize, ySize);
            temperature = std::make_shared<PWM::PWMDataStructure::flatStaggeredGrid<double>>(nX, nY, 0.5, 0.5, xSize, ySize);
            moisture = std::make_shared<PWM::PWMDataStructure::flatStaggeredGrid<double>>(nX, nY, 0.5, 0.5, xSize, ySize);
            ash = std::make_shared<PWM::PWMDataStructure::flatStaggeredGrid<double>>(nX, nY, 0.5, 0.5, xSize, ySize);
            cellTerrainType = std::make_shared<PWM::PWMDataStructure::flatStaggeredGrid<std::string>>(nX, nY, 0.5, 0.5, xSize, ySize);

            elevation->copy(other.getElevation());
            temperature->copy(other.getTemperature());
            moisture->copy(other.getMoisture());
            ash->copy(other.getMoisture());
            cellTerrainType->copy(other.getTerrainType());
        }

        template<typename T, typename TT, typename V, typename VV>
        inline terrain<T, TT, V, VV>& terrain<T, TT, V, VV>::operator=(terrain<T, TT, V, VV>& other){
            Planet = other.getPlanet();

            size_t width = other.getElevation().getWidth();
            V size = other.getElevation().getSize();
            elevation = std::make_shared<T>(T(width, size));
            temperature = std::make_shared<T>(T(width, size));
            moisture = std::make_shared<T>(T(width, size));
            ash = std::make_shared<T>(T(width, size));
            cellTerrainType = std::make_shared<TT>(TT(width, size));
            
            elevation->copy(other.getElevation());
            temperature->copy(other.getTemperature());
            moisture->copy(other.getMoisture());
            ash->copy(other.getAsh());
            cellTerrainType->copy(other.getTerrainType());

            return *this;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const T& terrain<T, TT, V, VV>::getElevation() const{
            return *elevation;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const T& terrain<T, TT, V, VV>::getTemperature() const{
            return *temperature;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const T& terrain<T, TT, V, VV>::getMoisture() const{
            return *moisture;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const T &terrain<T, TT, V, VV>::getAsh() const{
            return *ash;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const TT& terrain<T, TT, V, VV>::getTerrainType() const{
            return *cellTerrainType;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<planet> terrain<T, TT, V, VV>::getPlanet() const{
            return Planet;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V terrain<T, TT, V, VV>::getElevation(int index) const{
            return elevation->getData(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V terrain<T, TT, V, VV>::getTemperature(int index) const{
            return temperature->getData(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V terrain<T, TT, V, VV>::getMoisture(int index) const{
            return moisture->getData(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V terrain<T, TT, V, VV>::getAsh(int index) const{
            return ash->getData(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const VV terrain<T, TT, V, VV>::getTerrainType(int index) const{
            return cellTerrainType->getData(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::setTemperature(int index, V& newTemp){
            temperature->setData(index, newTemp);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::setMoisture(int index, V& newMoist){
            moisture->setData(index, newMoist);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::setAsh(int index, V& newAsh){
            ash->setData(index, newAsh);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::addAsh(int index, V &ashAdded){
            V newAsh = this->getAsh(index) + ashAdded;
            this->setAsh(index, newAsh);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::setTerrainType(int index, VV& newType){
            cellTerrainType->setData(index, newType);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::setElevation(int index, V& newElev){
            elevation->setData(index, newElev);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline bool terrain<T, TT, V, VV>::operator==(terrain<T, TT, V, VV>& other) const{
            if (*(getPlanet()) != *(other.getPlanet()))
                return false;
            if (getElevation() != other.getElevation())
                return false;
            if (getTemperature() != other.getTemperature())
                return false;
            if (getMoisture() != other.getMoisture())
                return false;
            if (getAsh() != other.getAsh())
                return false;
            if (getTerrainType() != other.getTerrainType())
                return false;
            return true;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline bool terrain<T, TT, V, VV>::operator!=(terrain<T, TT, V, VV>& other) const{
            return !(*this == other);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void terrain<T, TT, V, VV>::init(std::vector<std::pair<V, bool>>& terrainData){
            if (elevation->size() != terrainData.size()){
                std::cerr << "Error! Incorrect number of datapoints given for terrain initialisation!" << std::endl;
                std::cerr << "Elevation has size: " << elevation->size() << " compared to terrainData size: " << terrainData.size() << "!" << std::endl;
                return;
            }
            #pragma omp parallel for
                for (int i = 0; i < terrainData.size(); ++i){
                    V h = terrainData.at(i).first;
                    elevation->setData(i, h);
                    if (terrainData.at(i).second){
                        std::string terT = "WATER";
                        setTerrainType(i, terT);
                        V mst = 1000;
                        setMoisture(i, mst);
                    }
                    else{
                        std::string terT = "SOIL";
                        setTerrainType(i, terT);
                        V mst = 0;
                        setMoisture(i, mst);
                    }
                    V templess = (h <= 0) ? 273.15 : Planet->getAveTemp() - (h / 100);
                    setTemperature(i, templess);
                }
            #pragma omp barrier
        }

        template<>
        inline void terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::init(std::vector<std::pair<double, bool>>& terrainData){
            if (elevation->size() != terrainData.size()){
                std::cerr << "Error! Incorrect number of datapoints given for terrain initialisation!" << std::endl;
                std::cerr << "Elevation has size: " << elevation->size() << " compared to terrainData size: " << terrainData.size() << "!" << std::endl;
                return;
            }
            #pragma omp parallel for
                for (int i = 0; i < terrainData.size(); ++i){
                    double h = terrainData.at(i).first;
                    elevation->setData(i, h);
                    if (terrainData.at(i).second){
                        std::string terT = "WATER";
                        setTerrainType(i, terT);
                        double mst = 1000;
                        setMoisture(i, mst);
                    }
                    else{
                        std::string terT = "SOIL";
                        setTerrainType(i, terT);
                        double mst = 0;
                        setMoisture(i, mst);
                    }
                    double templess = (h <= 0) ? 273.15 : Planet->getAveTemp() - (h / 100);
                    setTemperature(i, templess);
                }
            #pragma omp barrier
        }

        template<>
        inline void terrain<PWM::PWMDataStructure::square2DArray<float>, PWM::PWMDataStructure::square2DArray<std::string>, float, std::string>::init(std::vector<std::pair<float, bool>>& terrainData){
            if (elevation->size() != terrainData.size()){
                std::cerr << "Error! Incorrect number of datapoints given for terrain initialisation!" << std::endl;
                std::cerr << "Elevation has size: " << elevation->size() << " compared to terrainData size: " << terrainData.size() << "!" << std::endl;
                return;
            }
            #pragma omp parallel for
                for (int i = 0; i < terrainData.size(); ++i){
                    float h = terrainData.at(i).first;
                    elevation->setData(i, h);
                    if (terrainData.at(i).second){
                        std::string terT = "WATER";
                        setTerrainType(i, terT);
                        float mst = 1000;
                        setMoisture(i, mst);
                    }
                    else{
                        std::string terT = "SOIL";
                        setTerrainType(i, terT);
                        float mst = 0;
                        setMoisture(i, mst);
                    }
                    float templess = (h <= 0) ? 273.15 : Planet->getAveTemp() - (h / 100);
                    setTemperature(i, templess);
                }
            #pragma omp barrier
        }
    }
}
#endif //PWM_MODEL_TERRAIN_H
