#ifndef PWM_MODEL_WORLD_H
#define PWM_MODEL_WORLD_H

#include "airLayer.h"
#include "convectionLayer.h"
#include <fstream>
#include <iostream>
#include <limits>
#include "planet.h"
#include "pressureWave.h"
#include "settings.h"
#include "square2DArray.h"
#include "sun.h"
#include "terrain.h"
#include "terrain_structure.hpp"
#include <vector>

namespace PWM{
    namespace Model{
        template<typename T, typename TT, typename V, typename VV> class world{
            private:
                std::shared_ptr<terrain<T, TT, V, VV>> Terrain;

                std::shared_ptr<planet> Planet;

                std::vector<std::shared_ptr<sun<V>>> suns;

                std::vector<std::shared_ptr<airLayer<T, V>>> airLayers;
                std::vector<std::shared_ptr<convectionLayer<T, V>>> convectionLayers;
                float threshold = 0.4f;

                V actualSize = 50000; //actual size of a side in metres
                size_t terrainWidth; //number of cells to a side

                size_t terWidth; //Number of cells in the x-direction
                size_t terHeight; //Number of cells in the y-direction
                double worldXLength; //Length of the world in the X direction, in metres.
                double worldYLength; //Length of the world in the Y direction, in metres.

                void sortLayers();
            public:
                world();
                world(std::string pFile, std::string terFile, const V actSize = 50000);
                world(std::string pFile, terrain_structure& terrStruct, const size_t width, const V actSize = 50000);
                world(std::string pFile, terrain_structure& terrStruct, const size_t width, const size_t height, const V xWrldSize, const V yWrldSize);
                world(PWM::Utils::settings& s);
                world(PWM::Utils::settings& s, terrain_structure& terrStruct);

                const std::shared_ptr<terrain<T, TT, V, VV>>& getTerrain() const;
                const std::shared_ptr<planet>& getPlanet() const;
                const std::vector<std::shared_ptr<airLayer<T, V>>>& getAirLayers() const;
                const std::shared_ptr<airLayer<T, V>>& getAirLayer(size_t index) const;
                const std::vector<std::shared_ptr<convectionLayer<T, V>>>& getConvectionLayers() const;
                const std::shared_ptr<convectionLayer<T, V>>& getConvectionLayer(size_t index) const;
                const std::vector<std::shared_ptr<sun<V>>>& getSuns() const;
                const std::shared_ptr<sun<V>>& getSun(const size_t index) const;
                const V getWorldSize() const;
                const V getCellSize() const;
                const V getWorldXSize() const;
                const V getWorldYSize() const;
                const size_t getWorldCellWidth() const;
                const size_t getWorldCellHeight() const;
                const std::pair<V, V> getWorldXYCoords(const int i, const int j);

                const std::pair<V, V> getLocalVel(const float x, const float y, const float z);
                const V getLocalVelX(const float x, const float y, const float z);
                const V getLocalVelY(const float x, const float y, const float z);
                const V getLocalTemp(const float x, const float y, const float z);
                const V getLocalMoist(const float x, const float y, const float z);
                const V getLocalClouds(const float x, const float y, const float z);
                const V getLocalConvection(const float x, const float y, const float z);

                const std::pair<V, V> getGlobalVel(const float x, const float y, const float z, const float radius);
                const V getGlobalVelX(const float x, const float y, const float z, const float radius);
                const V getGlobalVelY(const float x, const float y, const float z, const float radius);
                const V getGlobalTemp(const float x, const float y, const float z, const float radius);
                const V getGlobalMoist(const float x, const float y, const float z, const float radius);
                const V getGlobalClouds(const float x, const float y, const float z, const float radius);

                const int getLayerAtHeight(const V h);

                void init(std::string pFile, std::string terFile);
                void init(std::string pFile, terrain_structure terrStruct);
                void init(std::string pFile, terrain_structure terrStruct, size_t width);
                void init(std::string pFile, terrain_structure terrStruct, std::vector<V>& layerHeights);
                void init(std::string pFile, terrain_structure terrStruct, const size_t width, const size_t height, const V xWrldSize, const V yWrldSize);
                void init(std::string pFile, terrain_structure terrStruct, const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, std::vector<V>& layerHeights);
                void init(terrain_structure terrStruct);

                void addSun(std::string sFile);

                void assignAirLayerObstacles();

                bool operator==(const world<T, TT, V, VV>& rhs) const;
                int maxLayers = 10;
                float atmoTop = 10000.f;
                float maxThickness = 1000.f;
        };
        
        template<typename T, typename TT, typename V, typename VV>
        inline world<T, TT, V, VV>::world(){
            Terrain = nullptr;
            Planet = nullptr;
            airLayers = std::vector<std::shared_ptr<airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<convectionLayer<T, V>>>();
            suns = std::vector<std::shared_ptr<sun<V>>>();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline world<T, TT, V, VV>::world(std::string pFile, std::string terFile, const V actSize): actualSize(actSize){
            Terrain = nullptr;
            Planet = nullptr;
            airLayers = std::vector<std::shared_ptr<airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<convectionLayer<T, V>>>();
            suns = std::vector<std::shared_ptr<sun<V>>>();
            init(pFile, terFile);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline world<T, TT, V, VV>::world(std::string pFile, terrain_structure& terrStruct, const size_t width, const V actSize): terrainWidth(width), actualSize(actSize){
            Terrain = nullptr;
            Planet = nullptr;
            airLayers = std::vector<std::shared_ptr<airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<convectionLayer<T, V>>>();
            suns = std::vector<std::shared_ptr<sun<V>>>();
            init(pFile, terrStruct);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline world<T, TT, V, VV>::world(std::string pFile, terrain_structure &terrStruct, const size_t width, const size_t height, const V xWrldSize, const V yWrldSize) : terWidth(width), terHeight(height), worldXLength(xWrldSize), worldYLength(yWrldSize){
            Terrain = nullptr;
            Planet = nullptr;
            airLayers = std::vector<std::shared_ptr<airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<convectionLayer<T, V>>>();
            suns = std::vector<std::shared_ptr<sun<V>>>();
            init(pFile, terrStruct);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline world<T, TT, V, VV>::world(PWM::Utils::settings& s):    actualSize(s.actualSize),
                                                                threshold(s.wmThreshold),
                                                                atmoTop(s.atmoTopHeight),
                                                                maxLayers(s.maxLayers),
                                                                terrainWidth(s.terrainCellWidth){
            Terrain = nullptr;
            auto pl = s.P;
            Planet = std::make_shared<PWM::Model::planet>(pl);
            airLayers = std::vector<std::shared_ptr<airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<convectionLayer<T, V>>>();
            suns = std::vector<std::shared_ptr<sun<V>>>();
            for (auto x : s.suns){
                PWM::Model::sun<V> sn = x;
                suns.push_back(std::make_shared<PWM::Model::sun<V>>(sn));
            }
            terrain_structure terrStruct(s.terrainFile, s.gridWidth, s.actualSize);
            init(terrStruct, terrainWidth);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline world<T, TT, V, VV>::world(PWM::Utils::settings& s, terrain_structure& terrStruct):    actualSize(s.actualSize),
                                                                threshold(s.wmThreshold),
                                                                atmoTop(s.atmoTopHeight),
                                                                maxLayers(s.maxLayers),
                                                                terrainWidth(s.terrainCellWidth){
            Terrain = nullptr;
            auto pl = s.P;
            Planet = std::make_shared<PWM::Model::planet>(pl);
            airLayers = std::vector<std::shared_ptr<airLayer<T, V>>>();
            convectionLayers = std::vector<std::shared_ptr<convectionLayer<T, V>>>();
            suns = std::vector<std::shared_ptr<sun<V>>>();
            for (auto x : s.suns){
                PWM::Model::sun<V> sn = x;
                suns.push_back(std::make_shared<PWM::Model::sun<V>>(sn));
            }
            init(terrStruct, terrainWidth);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<terrain<T, TT, V, VV>>& world<T, TT, V, VV>::getTerrain() const{
            if (Terrain == nullptr)
                std::cerr << "Initialise the terrain first!" << std::endl;
            return Terrain;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<planet>& world<T, TT, V, VV>::getPlanet() const{
            if (Planet == nullptr)
                std::cerr << "Initialise the planet first!" << std::endl;
            return Planet;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<airLayer<T, V>>>& world<T, TT, V, VV>::getAirLayers() const{
            if (airLayers.size() == 0)
                std::cerr << "Initialise the world first!" << std::endl;
            return airLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<airLayer<T, V>>& world<T, TT, V, VV>::getAirLayer(size_t index) const{
            return airLayers.at(index);    
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<convectionLayer<T, V>>>& world<T, TT, V, VV>::getConvectionLayers() const{
            if (convectionLayers.size() == 0)
                std::cerr << "Initialise the world first!" << std::endl;
            return convectionLayers;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<convectionLayer<T, V>>& world<T, TT, V, VV>::getConvectionLayer(size_t index) const{
            return convectionLayers.at(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::vector<std::shared_ptr<sun<V>>>& world<T, TT, V, VV>::getSuns() const{
            return suns;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::shared_ptr<sun<V>>& world<T, TT, V, VV>::getSun(const size_t index) const{
            return suns.at(index);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getWorldSize() const{
            return actualSize;
        }

        template<>
        inline const double world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::getWorldSize() const{
            return worldXLength;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const size_t world<T, TT, V, VV>::getWorldCellWidth() const{
            return terrainWidth;
        }

        template<>
        inline const size_t world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::getWorldCellWidth() const{
            return terWidth;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const size_t world<T, TT, V, VV>::getWorldCellHeight() const{
            return terrainWidth;
        }

        template<>
        inline const size_t world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::getWorldCellHeight() const{
            return terHeight;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::pair<V, V> world<T, TT, V, VV>::getWorldXYCoords(const int i, const int j){
            int i1 = i - (getWorldCellWidth() / 2);
            int j1 = j - (getWorldCellWidth() / 2);
            V x = (i1 + 0.5) * getCellSize();
            V y = (j1 + 0.5) * getCellSize();
            return std::make_pair(x, y);
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getCellSize() const{
            return this->getWorldSize() / this->getWorldCellWidth();
        }

        template<>
        inline const double world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::getCellSize() const{
            return this->getWorldXSize() / this->getWorldCellWidth();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getWorldXSize() const{
            return worldXLength;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getWorldYSize() const{
            return worldYLength;
        }
        
        template<typename T, typename TT, typename V, typename VV>
        inline const std::pair<V, V> world<T, TT, V, VV>::getLocalVel(const float x, const float y, const float z){
            return std::make_pair(getLocalVelY(x, y, z), getLocalVelX(x, y, z));
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getLocalVelY(const float x, const float y, const float z){
            int firstTopIndex;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i)->getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            if (firstTopIndex == 0)
                return getAirLayer(firstTopIndex)->getVelocityTheta().getInterpolated(x, y);
            else if (firstTopIndex == airLayers.size())
                return getAirLayer(firstTopIndex - 1)->getVelocityTheta().getInterpolated(x, y);
            else {
                V factor = (getAirLayer(firstTopIndex)->getHeight() - z) / (getAirLayer(firstTopIndex)->getHeight() - getAirLayer(firstTopIndex - 1)->getHeight());
                return PWM::Utils::interpolate(getAirLayer(firstTopIndex)->getVelocityTheta().getInterpolated(x, y), getAirLayer(firstTopIndex - 1)->getVelocityTheta().getInterpolated(x, y), factor);
            }
        }
        
        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getLocalVelX(const float x, const float y, const float z){
            int firstTopIndex;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i)->getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            if (firstTopIndex == 0)
                return getAirLayer(firstTopIndex)->getVelocityPhi().getInterpolated(x, y);
            else if (firstTopIndex == airLayers.size())
                return getAirLayer(firstTopIndex - 1)->getVelocityPhi().getInterpolated(x, y);
            else {
                V factor = (getAirLayer(firstTopIndex)->getHeight() - z) / (getAirLayer(firstTopIndex)->getHeight() - getAirLayer(firstTopIndex - 1)->getHeight());
                return PWM::Utils::interpolate(getAirLayer(firstTopIndex)->getVelocityPhi().getInterpolated(x, y), getAirLayer(firstTopIndex - 1)->getVelocityPhi().getInterpolated(x, y), factor);
            }
        }
        
        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getLocalTemp(const float x, const float y, const float z){
            int firstTopIndex;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i)->getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            if (firstTopIndex == 0)
                return getAirLayer(firstTopIndex)->getTemperature().getInterpolated(x, y);
            else if (firstTopIndex == airLayers.size())
                return getAirLayer(firstTopIndex - 1)->getTemperature().getInterpolated(x, y);
            else {
                V factor = (getAirLayer(firstTopIndex)->getHeight() - z) / (getAirLayer(firstTopIndex)->getHeight() - getAirLayer(firstTopIndex - 1)->getHeight());
                return PWM::Utils::interpolate(getAirLayer(firstTopIndex)->getTemperature().getInterpolated(x, y), getAirLayer(firstTopIndex - 1)->getTemperature().getInterpolated(x, y), factor);
            }
        }
        
        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getLocalMoist(const float x, const float y, const float z){
            int firstTopIndex;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i)->getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            if (firstTopIndex == 0)
                return getAirLayer(firstTopIndex)->getMoisture().getInterpolated(x, y);
            else if (firstTopIndex == airLayers.size())
                return getAirLayer(firstTopIndex - 1)->getMoisture().getInterpolated(x, y);
            else {
                V factor = (getAirLayer(firstTopIndex)->getHeight() - z) / (getAirLayer(firstTopIndex)->getHeight() - getAirLayer(firstTopIndex - 1)->getHeight());
                return PWM::Utils::interpolate(getAirLayer(firstTopIndex)->getMoisture().getInterpolated(x, y), getAirLayer(firstTopIndex - 1)->getMoisture().getInterpolated(x, y), factor);
            }
        }
        
        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getLocalClouds(const float x, const float y, const float z){
            int firstTopIndex;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i)->getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            if (firstTopIndex == 0)
                return getAirLayer(firstTopIndex)->getCondensedWater().getInterpolated(x, y);
            else if (firstTopIndex == airLayers.size())
                return getAirLayer(firstTopIndex - 1)->getCondensedWater().getInterpolated(x, y);
            else {
                V factor = (getAirLayer(firstTopIndex)->getHeight() - z) / (getAirLayer(firstTopIndex)->getHeight() - getAirLayer(firstTopIndex - 1)->getHeight());
                return PWM::Utils::interpolate(getAirLayer(firstTopIndex)->getCondensedWater().getInterpolated(x, y), getAirLayer(firstTopIndex - 1)->getCondensedWater().getInterpolated(x, y), factor);
            }
        }
    
        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getLocalConvection(const float x, const float y, const float z){
            int belowIndex = 0;
            for(int i = convectionLayers.size() - 1; i >= 0; i--){
                if(z > getConvectionLayer(i)->getLayerHeight()){
                    belowIndex = i;
                }
            }
            if (belowIndex == 0 || belowIndex == convectionLayers.size() - 1){
                return getConvectionLayer(belowIndex)->getVerticalVelocities().getInterpolated(x,y);
            }
            V factor = (z - getConvectionLayer(belowIndex)->getLayerHeight()) / (getConvectionLayer(belowIndex + 1)->getLayerHeight() - getConvectionLayer(belowIndex)->getLayerHeight());
            return (getConvectionLayer(belowIndex)->getVerticalVelocities().getInterpolated(x,y), getConvectionLayer(belowIndex + 1)->getVerticalVelocities().getInterpolated(x,y), factor); // JG This line is broken Cilliers to fix
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const std::pair<V, V> world<T, TT, V, VV>::getGlobalVel(const float x, const float y, const float z, const float radius){
            return std::make_pair(getGlobalVelY(x, y, z, radius), getGlobalVelX(x, y, z, radius));
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getGlobalVelY(const float x, const float y, const float z, const float radius){
            int firstTopIndex = 0;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i).getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            V topSum = 0, botSum = 0;
            int cellCount = 0;
            for (int i = floor(x - radius); i < (x + radius); ++i){
                for (int j = floor(y - radius); j < (y + radius); ++j){
                    ++cellCount;
                    if (firstTopIndex == 0)
                        topSum += getAirLayer(firstTopIndex).getVelocityTheta(i, j);
                    else if (firstTopIndex == airLayers.size()){
                        botSum += getAirLayer(firstTopIndex - 1).getVelocityTheta(i, j);
                    }
                    else{
                        botSum += getAirLayer(firstTopIndex).getVelocityTheta(i, j);
                        topSum += getAirLayer(firstTopIndex + 1).getVelocityTheta(i, j);
                    }
                }
            }

            if (firstTopIndex == 0){
                return topSum / cellCount;
            }
            else if (firstTopIndex == airLayers.size()){
                return botSum / cellCount;
            }
            else{
                V factor = (getAirLayer(firstTopIndex) - z) / (getAirLayer(firstTopIndex).getHeight() - getAirLayer(firstTopIndex - 1).getHeight());
                return PWM::Utils::interpolate(topSum / cellCount, botSum / cellCount);
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getGlobalVelX(const float x, const float y, const float z, const float radius){
            int firstTopIndex = 0;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i).getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            V topSum = 0, botSum = 0;
            int cellCount = 0;
            for (int i = floor(x - radius); i < (x + radius); ++i){
                for (int j = floor(y - radius); j < (y + radius); ++j){
                    ++cellCount;
                    if (firstTopIndex == 0)
                        topSum += getAirLayer(firstTopIndex).getVelocityPhi(i, j);
                    else if (firstTopIndex == airLayers.size()){
                        botSum += getAirLayer(firstTopIndex - 1).getVelocityPhi(i, j);
                    }
                    else{
                        botSum += getAirLayer(firstTopIndex).getVelocityPhi(i, j);
                        topSum += getAirLayer(firstTopIndex + 1).getVelocityPhi(i, j);
                    }
                }
            }

            if (firstTopIndex == 0){
                return topSum / cellCount;
            }
            else if (firstTopIndex == airLayers.size()){
                return botSum / cellCount;
            }
            else{
                V factor = (getAirLayer(firstTopIndex) - z) / (getAirLayer(firstTopIndex).getHeight() - getAirLayer(firstTopIndex - 1).getHeight());
                return PWM::Utils::interpolate(topSum / cellCount, botSum / cellCount);
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getGlobalTemp(const float x, const float y, const float z, const float radius){
            int firstTopIndex = 0;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i).getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            V topSum = 0, botSum = 0;
            int cellCount = 0;
            for (int i = floor(x - radius); i < (x + radius); ++i){
                for (int j = floor(y - radius); j < (y + radius); ++j){
                    ++cellCount;
                    if (firstTopIndex == 0)
                        topSum += getAirLayer(firstTopIndex).getTemperature(i, j);
                    else if (firstTopIndex == airLayers.size()){
                        botSum += getAirLayer(firstTopIndex - 1).getTemperature(i, j);
                    }
                    else{
                        botSum += getAirLayer(firstTopIndex).getTemperature(i, j);
                        topSum += getAirLayer(firstTopIndex + 1).getTemperature(i, j);
                    }
                }
            }

            if (firstTopIndex == 0){
                return topSum / cellCount;
            }
            else if (firstTopIndex == airLayers.size()){
                return botSum / cellCount;
            }
            else{
                V factor = (getAirLayer(firstTopIndex) - z) / (getAirLayer(firstTopIndex).getHeight() - getAirLayer(firstTopIndex - 1).getHeight());
                return PWM::Utils::interpolate(topSum / cellCount, botSum / cellCount);
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getGlobalMoist(const float x, const float y, const float z, const float radius){
            int firstTopIndex = 0;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i).getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            V topSum = 0, botSum = 0;
            int cellCount = 0;
            for (int i = floor(x - radius); i < (x + radius); ++i){
                for (int j = floor(y - radius); j < (y + radius); ++j){
                    ++cellCount;
                    if (firstTopIndex == 0)
                        topSum += getAirLayer(firstTopIndex).getMoisture(i, j);
                    else if (firstTopIndex == airLayers.size()){
                        botSum += getAirLayer(firstTopIndex - 1).getMoisture(i, j);
                    }
                    else{
                        botSum += getAirLayer(firstTopIndex).getMoisture(i, j);
                        topSum += getAirLayer(firstTopIndex + 1).getMoisture(i, j);
                    }
                }
            }

            if (firstTopIndex == 0){
                return topSum / cellCount;
            }
            else if (firstTopIndex == airLayers.size()){
                return botSum / cellCount;
            }
            else{
                V factor = (getAirLayer(firstTopIndex) - z) / (getAirLayer(firstTopIndex).getHeight() - getAirLayer(firstTopIndex - 1).getHeight());
                return PWM::Utils::interpolate(topSum / cellCount, botSum / cellCount);
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const V world<T, TT, V, VV>::getGlobalClouds(const float x, const float y, const float z, const float radius){
            int firstTopIndex = 0;
            for(int i = 0; i < airLayers.size(); ++i){
                if(z < getAirLayer(i).getHeight()){
                    firstTopIndex = i;
                    break;
                }
                else if (i == airLayers.size() - 1){
                    firstTopIndex = i + 1;
                }
            }
            V topSum = 0, botSum = 0;
            int cellCount = 0;
            for (int i = floor(x - radius); i < (x + radius); ++i){
                for (int j = floor(y - radius); j < (y + radius); ++j){
                    ++cellCount;
                    if (firstTopIndex == 0)
                        topSum += getAirLayer(firstTopIndex).getCondensedWater(i, j);
                    else if (firstTopIndex == airLayers.size()){
                        botSum += getAirLayer(firstTopIndex - 1).getCondensedWater(i, j);
                    }
                    else{
                        botSum += getAirLayer(firstTopIndex).getCondensedWater(i, j);
                        topSum += getAirLayer(firstTopIndex + 1).getCondensedWater(i, j);
                    }
                }
            }

            if (firstTopIndex == 0){
                return topSum / cellCount;
            }
            else if (firstTopIndex == airLayers.size()){
                return botSum / cellCount;
            }
            else{
                V factor = (getAirLayer(firstTopIndex) - z) / (getAirLayer(firstTopIndex).getHeight() - getAirLayer(firstTopIndex - 1).getHeight());
                return PWM::Utils::interpolate(topSum / cellCount, botSum / cellCount);
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        inline const int world<T, TT, V, VV>::getLayerAtHeight(const V h){
            for (int i = 0; i < airLayers.size(); ++i){
                if (airLayers[i]->getLayerTop() > h)
                    return i;
            }
            return airLayers.size() - 1;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(std::string pFile, std::string terFile){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain file \"" << terFile << "\" halted." << std::endl;
        }

        /*template<>
        void world<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::init(std::string pFile, std::string terFile){
            Planet = std::make_shared<planet>(planet(pFile));
            auto fil = std::fstream(terFile);
            if (!(fil.good())){
                std::cerr << "Error with opening the terrain file input!" << std::endl;
                return;
            }
            double minLandElev = std::numeric_limits<double>::max();
            int t, p;
            std::vector<std::pair<double, bool>> datapoints;
            try{
                double min, max;
                fil >> t;
                fil >> p;
                fil >> min;
                fil >> max;
                terrain<PWM::PWMDataStructure::SSGDataStructure<double>, PWM::PWMDataStructure::SSGDataStructure<std::string>, double, std::string> x = terrain<PWM::PWMDataStructure::SSGDataStructure<double>, PWM::PWMDataStructure::SSGDataStructure<std::string>, double, std::string>(Planet, t, p);
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::SSGDataStructure<double>, PWM::PWMDataStructure::SSGDataStructure<std::string>, double, std::string>>(x);
                double elev;
                bool ocean;
                int counter = 0;
                while (counter < t * p){
                    fil >> elev;
                    fil >> ocean;
                    if (!ocean && elev < minLandElev)
                        minLandElev = elev;
                    datapoints.push_back(std::make_pair(elev, ocean));
                    ++counter;
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }
            auto atGround = [](const float threshold, const double elevation, const std::vector<std::pair<double, bool>>& data)->bool{
                int count = 0;
                float thres = threshold * data.size();
                for (auto x : data){
                    if (x.first >= elevation){
                        ++count;
                        if (count >= thres)
                            return true;
                    }
                }
                return false;
            };

            double airLayerThickness = atmoTop / maxLayers;
            double alt = minLandElev + (airLayerThickness / 2);
            while (alt < atmoTop && airLayers.size() < maxLayers){
                auto x = std::make_shared<airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(airLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, alt, airLayerThickness, getWorldCellWidth(), getWorldSize()));
                x->randomInit();
                airLayers.push_back(x);
                alt += airLayerThickness;
            }
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), t, p));
                convectionLayers.push_back(x);
            }
        }*/

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(std::string pFile, terrain_structure terrStruct){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain struct halted." << std::endl;
        }

        template<>
        inline void world<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::init(std::string pFile, terrain_structure terrStruct){
            Planet = std::make_shared<planet>(planet(pFile));
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
            if (getWorldCellWidth() != terrStruct.field_size)
                terrainWidth = std::pow(2, std::ceil(std::log2(terrStruct.field_size)));
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
            }
            try{
//                double min, max;
                terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>(Planet, getWorldCellWidth(), getWorldSize());
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
//                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < getWorldCellWidth(); ++i){
                    for (int j = 0; j < getWorldCellWidth(); ++j){
                        float x_i = i * terrStruct.field_size / getWorldCellWidth();
                        float y_j = j * terrStruct.field_size / getWorldCellWidth();
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            double airLayerThickness = atmoTop / maxLayers;
            double alt = minLandElev + (airLayerThickness / 2);
            while (alt < atmoTop && airLayers.size() < maxLayers){
                auto x = std::make_shared<airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(airLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, alt, airLayerThickness, getWorldCellWidth(), getWorldSize()));
                x->randomInit();
                airLayers.push_back(x);
                alt += airLayerThickness;
            }
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), getWorldCellWidth(), getWorldSize()));
                convectionLayers.push_back(x);
            }
            assignAirLayerObstacles();
        }

        template<>
        inline void world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::init(std::string pFile, terrain_structure terrStruct){
            Planet = std::make_shared<planet>(planet(pFile));
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
            if (getWorldCellWidth() != terrStruct.field_size){
                terrainWidth = terrStruct.field_size;
                terWidth = terrStruct.field_size;
                terHeight = terrStruct.field_size;
            }
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
                worldXLength = interval_size;
                worldYLength = interval_size;
            }
            try{
//                double min, max;
                terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>(Planet, terWidth, terHeight, worldXLength, worldYLength);
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
//                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < terHeight; ++i){
                    for (int j = 0; j < terWidth; ++j){
                        float x_i = i * terrStruct.field_size / worldXLength;
                        float y_j = j * terrStruct.field_size / worldYLength;
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            double airLayerThickness = atmoTop / maxLayers;
            double alt = minLandElev + (airLayerThickness / 2);
            while (alt < atmoTop && airLayers.size() < maxLayers){
                auto x = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, alt, airLayerThickness, terWidth, terHeight, worldXLength, worldYLength));
                x->randomInit();
                airLayers.push_back(x);
                alt += airLayerThickness;
            }
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), terWidth, terHeight, worldXLength, worldYLength));
                convectionLayers.push_back(x);
            }
            assignAirLayerObstacles();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(std::string pFile, terrain_structure terrStruct, size_t width){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain struct halted." << std::endl;
        }

        template<>
        inline void world<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::init(std::string pFile, terrain_structure terrStruct, size_t width){
            Planet = std::make_shared<planet>(planet(pFile));
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
            terrainWidth = std::pow(2, std::ceil(std::log2(width)));
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
            }
            try{
//                double min, max;
                terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>(Planet, getWorldCellWidth(), getWorldSize());
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
//                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < getWorldCellWidth(); ++i){
                    for (int j = 0; j < getWorldCellWidth(); ++j){
                        float x_i = i * terrStruct.field_size / getWorldCellWidth();
                        float y_j = j * terrStruct.field_size / getWorldCellWidth();
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            double airLayerThickness = atmoTop / maxLayers;
            double alt = minLandElev + (airLayerThickness / 2);
            while (alt < atmoTop && airLayers.size() < maxLayers){
                auto x = std::make_shared<airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(airLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, alt, airLayerThickness, getWorldCellWidth(), getWorldSize()));
                x->randomInit();
                airLayers.push_back(x);
                alt += airLayerThickness;
            }
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), getWorldCellWidth(), getWorldSize()));
                convectionLayers.push_back(x);
            }
            assignAirLayerObstacles();
        }


        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(std::string pFile, terrain_structure terrStruct, const size_t width, const size_t height, const V xWrldSize, const V yWrldSize){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain struct halted." << std::endl;
        }

        template<>
        inline void world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::init(std::string pFile, terrain_structure terrStruct, const size_t width, const size_t height, const double xWrldSize, const double yWrldSize){
            terWidth = width;
            terHeight = height;
            worldXLength = xWrldSize;
            worldYLength = yWrldSize;
            Planet = std::make_shared<planet>(planet(pFile));
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
//            terrainWidth = std::pow(2, std::ceil(std::log2(width)));
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
            }
            try{
//                double min, max;
                terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>(Planet, width, height, xWrldSize, yWrldSize);
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
//                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < height; ++i){
                    for (int j = 0; j < width; ++j){
                        float x_i = i * terrStruct.field_size / height;
                        float y_j = j * terrStruct.field_size / width;
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            double airLayerThickness = atmoTop / maxLayers;
            double alt = minLandElev + (airLayerThickness / 2);
            while (alt < atmoTop && airLayers.size() < maxLayers){
                auto x = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, alt, airLayerThickness, width, height, xWrldSize, yWrldSize));
                x->randomInit();
                airLayers.push_back(x);
                alt += airLayerThickness;
            }
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), width, height, xWrldSize, yWrldSize));
                convectionLayers.push_back(x);
            }
            assignAirLayerObstacles();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(std::string pFile, terrain_structure terrStruct, const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, std::vector<V>& layerHeights){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain struct halted." << std::endl;
        }

        template<>
        inline void world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::init(std::string pFile, terrain_structure terrStruct, const size_t width, const size_t height, const double xWrldSize, const double yWrldSize, std::vector<double>& layerHeights){
            terWidth = width;
            terHeight = height;
            worldXLength = xWrldSize;
            worldYLength = yWrldSize;
            Planet = std::make_shared<planet>(planet(pFile));
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
//            terrainWidth = std::pow(2, std::ceil(std::log2(width)));
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
            }

            try{
//                double min, max;
                terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>(Planet, width, height, xWrldSize, yWrldSize);
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
//                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < height; ++i){
                    for (int j = 0; j < width; ++j){
                        float x_i = i * terrStruct.field_size / height;
                        float y_j = j * terrStruct.field_size / width;
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            if (layerHeights.size() > 0){
                //Initialise airlayers based on list of heights
                std::vector<std::pair<double, double>> res = std::vector<std::pair<double, double>>();
                res = PWM::Utils::allocateLayers<double>(layerHeights, res, maxThickness);

                double lHeight, thick;
                for (auto i : res){
                    thick = i.second - i.first;
                    lHeight = std::round((i.first + i.second) / 2);
                    auto y = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, lHeight, thick, width, height, xWrldSize, yWrldSize));
                    y->randomInit();
                    airLayers.push_back(y);
                }
                std::sort(airLayers.begin(), airLayers.end());
                //Fill in above the given layers (if there is space left)
                double remAtmo = atmoTop - airLayers.back()->getLayerTop();
                if (remAtmo >= 1 && airLayers.size() < maxLayers){
                    double thickness = remAtmo / (maxLayers - airLayers.size());
                    while (airLayers.size() < maxLayers){
                        double layerH = std::round(airLayers.back()->getLayerTop() + thickness / 2.0);
                        auto z = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, layerH, thick, width, height, xWrldSize, yWrldSize));
                        z->randomInit();
                        airLayers.push_back(z);
                    }
                }
            }
            else{
                double airLayerThickness = atmoTop / maxLayers;
                double alt = minLandElev + (airLayerThickness / 2);
                while (alt < atmoTop && airLayers.size() < maxLayers){
                    auto x = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, alt, airLayerThickness, width, height, xWrldSize, yWrldSize));
                    x->randomInit();
                    airLayers.push_back(x);
                    alt += airLayerThickness;
                }
            }
            //Add convection layers
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), width, height, xWrldSize, yWrldSize));
                convectionLayers.push_back(x);
            }
            //Paint the terrain obstacles into the air cells
            assignAirLayerObstacles();

        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(std::string pFile, terrain_structure terrStruct, std::vector<V> &layerHeights){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain struct halted." << std::endl;
        }

        template<>
        inline void world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::init(std::string pFile, terrain_structure terrStruct, std::vector<double>& layerHeights){
            terWidth = terrStruct.field_size;
            terHeight = terrStruct.field_size;
            worldXLength = terrStruct.max_xyz - terrStruct.min_xyz;
            worldYLength = terrStruct.max_xyz - terrStruct.min_xyz;
            Planet = std::make_shared<planet>(planet(pFile));
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
            //            terrainWidth = std::pow(2, std::ceil(std::log2(width)));
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
            }

            try{
                //                double min, max;
                terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>(Planet, terWidth, terHeight, worldXLength, worldYLength);
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
                //                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < terHeight; ++i){
                    for (int j = 0; j < terWidth; ++j){
                        float x_i = i * terrStruct.field_size / terHeight;
                        float y_j = j * terrStruct.field_size / terWidth;
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            if (layerHeights.size() > 0){
                //Initialise airlayers based on list of heights
                std::vector<std::pair<double, double>> res = std::vector<std::pair<double, double>>();
                res = PWM::Utils::allocateLayers<double>(layerHeights, res, maxThickness);

                double lHeight, thick;
                for (auto i : res){
                    thick = i.second - i.first;
                    lHeight = std::round((i.first + i.second) / 2);
                    auto y = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, lHeight, thick, terWidth, terHeight, worldXLength, worldYLength));
                    y->randomInit();
                    airLayers.push_back(y);
                }
                std::sort(airLayers.begin(), airLayers.end());
                //Fill in above the given layers (if there is space left)
                double remAtmo = atmoTop - airLayers.back()->getLayerTop();
                if (remAtmo >= 1 && airLayers.size() < maxLayers){
                    double thickness = remAtmo / (maxLayers - airLayers.size());
                    while (airLayers.size() < maxLayers){
                        double layerH = std::round(airLayers.back()->getLayerTop() + thickness / 2.0);
                        auto z = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, layerH, thick, terWidth, terHeight, worldXLength, worldYLength));
                        z->randomInit();
                        airLayers.push_back(z);
                    }
                }
            }
            else{
                double airLayerThickness = atmoTop / maxLayers;
                double alt = minLandElev + (airLayerThickness / 2);
                while (alt < atmoTop && airLayers.size() < maxLayers){
                    auto x = std::make_shared<airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, alt, airLayerThickness, terWidth, terHeight, worldXLength, worldYLength));
                    x->randomInit();
                    airLayers.push_back(x);
                    alt += airLayerThickness;
                }
            }
            //Add convection layers
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>(convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), terWidth, terHeight, worldXLength, worldYLength));
                convectionLayers.push_back(x);
            }
            //Paint the terrain obstacles into the air cells
            assignAirLayerObstacles();

        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::init(terrain_structure terrStruct){
            std::cerr << "Error! No data structure specified, initialisation of world with terrain struct halted." << std::endl;
        }

        template<>
        inline void world<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>::init(terrain_structure terrStruct){
            double minLandElev = std::numeric_limits<double>::max();
            std::vector<std::pair<double, bool>> datapoints;
            float interval_size = terrStruct.max_xyz - terrStruct.min_xyz;
            if (getWorldCellWidth() != terrStruct.field_size)
                terrainWidth = std::pow(2, std::ceil(std::log2(terrStruct.field_size)));
            if (getWorldSize() != interval_size){
                actualSize = interval_size;
            }
            try{
//                double min, max;
                terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string> temp = terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>(Planet, getWorldCellWidth(), getWorldSize());
                Terrain = std::make_shared<terrain<PWM::PWMDataStructure::square2DArray<double>, PWM::PWMDataStructure::square2DArray<std::string>, double, std::string>>(temp);
                double elev;
                bool ocean = false;
//                int x, y;
                auto mM = std::make_pair<double, double>(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
                for (int i = 0; i < terrStruct.height_field.size(); ++i){
                    for (int j = 0; j < terrStruct.height_field[0].size(); ++j){
                        if (terrStruct.height_field[i][j] < mM.first)
                            mM.first = terrStruct.height_field[i][j];
                    }
                }
                for (int i = 0; i < getWorldCellWidth(); ++i){
                    for (int j = 0; j < getWorldCellWidth(); ++j){
                        float x_i = i * terrStruct.field_size / getWorldCellWidth();
                        float y_j = j * terrStruct.field_size / getWorldCellWidth();
                        float x, y, s, t;
                        s = std::modf(x_i, &x);
                        t = std::modf(y_j, &y);
                        if (((int) x == terrStruct.field_size - 1) && ((int) y == terrStruct.field_size - 1)){
                            elev = terrStruct.height_field[x][y] - mM.first;
                        }
                        else if ((int) x == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_xy1, t);
                        }
                        else if ((int) y == terrStruct.field_size - 1){
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, s);
                        }
                        else{
                            double a_xy = terrStruct.height_field[(int) x][(int) y] - mM.first;
                            double a_x1y = terrStruct.height_field[(int) x + 1][(int) y] - mM.first;
                            double a_xy1 = terrStruct.height_field[(int) x][(int) y + 1] - mM.first;
                            double a_x1y1 = terrStruct.height_field[(int) x + 1][(int) y + 1] - mM.first;
                            elev = PWM::Utils::interpolate(a_xy, a_x1y, a_xy1, a_x1y1, s, t);
                        }
                        if (!ocean && elev < minLandElev)
                            minLandElev = elev;
                        datapoints.push_back(std::make_pair(elev, ocean));
                    }
                }
                Terrain->init(datapoints);
            }
            catch (std::exception& e){
                std::cerr << "Error when trying to parse the terrain file!" << std::endl;
                std::cerr << e.what() << std::endl;
                return;
            }

            double airLayerThickness = atmoTop / maxLayers;
            double alt = minLandElev + (airLayerThickness / 2);
            while (alt < atmoTop && airLayers.size() < maxLayers){
                auto x = std::make_shared<airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(airLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, alt, airLayerThickness, getWorldCellWidth(), getWorldSize()));
                x->randomInit();
                airLayers.push_back(x);
                alt += airLayerThickness;
            }
            sortLayers();
            for (int i = 0; i < airLayers.size() - 1; ++i){
                auto x = std::make_shared<convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>>(convectionLayer<PWM::PWMDataStructure::square2DArray<double>, double>(Planet, airLayers.at(i)->getHeight(), airLayers.at(i + 1)->getHeight(), getWorldCellWidth(), getWorldSize()));
                convectionLayers.push_back(x);
            }
            assignAirLayerObstacles();
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::sortLayers(){
            std::sort(airLayers.begin(), airLayers.end(), [](std::shared_ptr<PWM::Model::airLayer<T, V>>& l, std::shared_ptr<PWM::Model::airLayer<T, V>>& r){
                return l->getHeight() < r->getHeight();
            });
            std::sort(convectionLayers.begin(), convectionLayers.end(), [](std::shared_ptr<PWM::Model::convectionLayer<T, V>>& l, std::shared_ptr<PWM::Model::convectionLayer<T, V>>& r){
                return l->getLayerHeight() < r->getLayerHeight();
            });
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::assignAirLayerObstacles(){
            sortLayers();
            V zero = 0, one = 1;
            #pragma omp parallel for
                for (int i = 0; i < Terrain->getElevation().size(); ++i){
                    auto x = Terrain->getElevation(i);
                    for (int j = 0; j < airLayers.size(); ++j){
                        if (airLayers[j]->getLayerTop() >= x)
                            airLayers[j]->getObstacles().setData(i, zero);
                        else
                            airLayers[j]->getObstacles().setData(i, one);
                    }
                }
            #pragma omp barrier
        }

        template<typename T, typename TT, typename V, typename VV>
        inline bool world<T, TT, V, VV>::operator==(const world<T, TT, V, VV>& rhs) const{
            if (*getPlanet() != *(rhs.getPlanet()))
                return false;
            if (*getTerrain() != *(rhs.getTerrain()))
                return false;
            if (getAirLayers().size() != rhs.getAirLayers().size())
                return false;
            if (getConvectionLayers().size() != rhs.getConvectionLayers().size())
                return false;
            return true;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline void world<T, TT, V, VV>::addSun(const std::string sFile){
            auto s = std::make_shared<sun<V>>(sun<V>(sFile));
            suns.push_back(s);
        }
    }
}

#endif //PWM_MODEL_WORLD_H
