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
#ifndef PWM_UTILS_VDBEXPORTER_H
#define PWM_UTILS_VDBEXPORTER_H

#include <memory>
#include "openvdb/openvdb.h"
#include "world.h"

namespace PWM {
    namespace Utils {
        template <typename T, typename TT, typename V, typename VV> class vdbExporter{
            public:
                vdbExporter();
                int load(std::string filename, std::shared_ptr<PWM::Model::world<T, TT, V, VV>>);
                int write(const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputDirectory, int step, double sim_time, V cloudThres);
                int write(const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& airLayers, const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& convectionLayers, std::string outputDirectory, int step, double sim_time, V cloudThres);
        };

        template<typename T, typename TT, typename V, typename VV>
        vdbExporter<T, TT, V, VV>::vdbExporter(){
            openvdb::initialize();
        }

        template<typename T, typename TT, typename V, typename VV>
        int vdbExporter<T, TT, V, VV>::write(const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputDirectory, int step, double sim_time, V cloudThres){
            openvdb::GridPtrVec grids;
            openvdb::FloatGrid::Ptr cloudGrid = openvdb::FloatGrid::create();
            openvdb::Int32Grid::Ptr cloudTypeGrid = openvdb::Int32Grid::create();
            openvdb::FloatGrid::Ptr particulateGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr rainfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr ashfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velXGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velYGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velZGrid = openvdb::FloatGrid::create();
            cloudGrid->setName("cloudDensityGrid");
            grids.push_back(cloudGrid);
            cloudTypeGrid->setName("cloudTypeGrid");
            grids.push_back(cloudTypeGrid);
            particulateGrid->setName("particulateDensityGrid");
            grids.push_back(particulateGrid);
            rainfallGrid->setName("rainfallDensityGrid");
            grids.push_back(rainfallGrid);
            ashfallGrid->setName("ashfallDensityGrid");
            grids.push_back(ashfallGrid);
            velXGrid->setName("velXGrid");
            grids.push_back(velXGrid);
            velYGrid->setName("velYGrid");
            grids.push_back(velYGrid);
            velZGrid->setName("velZGrid");
            grids.push_back(velZGrid);

            auto cloud = cloudGrid->getAccessor();
            auto cType = cloudTypeGrid->getAccessor();
            auto parts = particulateGrid->getAccessor();
            auto rain = rainfallGrid->getAccessor();
            auto ash = ashfallGrid->getAccessor();
            auto velX = velXGrid->getAccessor();
            auto velY = velYGrid->getAccessor();
            auto velZ = velZGrid->getAccessor();
            int cellsPerSide = wld->getWorldCellWidth();
            auto cellWidth = wld->getCellSize();
            auto layers = wld->getAirLayers();
            auto convectLayers = wld->getConvectionLayers();
            V h, cH;
            for (int k = 0; k < layers.size(); ++k){
                h = layers[k]->getHeight();
                cloudGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                cloudGrid->insertMeta("cloudThreshold", openvdb::DoubleMetadata(cloudThres));
                cloudTypeGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                particulateGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velXGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velYGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                if (k < convectLayers.size()){
                    cH = convectLayers[k]->getLayerHeight();
                    velZGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    rainfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    ashfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                }
                for (int i = 0; i < cellsPerSide; ++i){
                    for (int j = 0; j < cellsPerSide; ++j){
//                        auto xy = wld->getWorldXYCoords(i, j);
                        openvdb::Coord ZYX(k, cellsPerSide - 1 - i, j);
                        cloud.setValue(ZYX, layers[k]->getCondensedWater(i, j));
                        cType.setValue(ZYX, layers[k]->getClouds(i, j));
                        parts.setValue(ZYX, layers[k]->getParticulates(i, j));
                        velX.setValue(ZYX, layers[k]->getVelocityPhi(i, j));
                        velY.setValue(ZYX, -(layers[k]->getVelocityTheta(i, j)));
                        if (k < convectLayers.size()){
                            rain.setValue(ZYX, convectLayers[k]->getRainfall(i, j));
                            velZ.setValue(ZYX, convectLayers[k]->getVerticalVelocity(i, j));
                            ash.setValue(ZYX, convectLayers[k]->getAshfall(i, j));
                        }
                    }
                }
            }
            for (auto x : grids){
                x->insertMeta("secondsSinceStart", openvdb::DoubleMetadata(sim_time));
                x->insertMeta("cellsPerSide", openvdb::Int32Metadata(cellsPerSide));
                x->insertMeta("cellSize", openvdb::DoubleMetadata(cellWidth));
            }

            
            std::stringstream fG;
            fG << outputDirectory << "/Weather_At_Frame_" << step << ".vdb";
            try {
                std::cerr << "Writing weather data to " << fG.str() << "." << std::endl;
                openvdb::io::File f(fG.str());
//                std::cerr << "File " << fG.str() << " has compression " << f.compression() << std::endl;
//                f.setCompression(openvdb::io::COMPRESS_ZIP | openvdb::io::COMPRESS_BLOSC);
                f.write(grids);
                f.close();
//                double i = 0;
//                for (auto x : grids){
//                    i += x->memUsage();
//                }

                for (auto x : grids)
                    x->clear();

//                double j = 0;
//                for (auto x : grids){
//                    j += x->memUsage();
//                }
//                std::cerr << "Mem usage before clear: " << i << std::endl;
//                std::cerr << "Mem usage after clear: " << j << std::endl;
                return 0;
            } catch (std::exception& e){
                std::cerr << "\033[1;31mError! Can't write weather data to file " << fG.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mWeather data not written.\033[0m" << std::endl;
                for (auto x : grids){
                    x->clear();
                }
                return -1;
            }
        }

        template<>
        int vdbExporter<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::write(const std::shared_ptr<PWM::Model::world<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>>& wld, std::string outputDirectory, int step, double sim_time, double cloudThres){
            openvdb::GridPtrVec grids;
            openvdb::FloatGrid::Ptr cloudGrid = openvdb::FloatGrid::create();
            openvdb::Int32Grid::Ptr cloudTypeGrid = openvdb::Int32Grid::create();
            openvdb::FloatGrid::Ptr particulateGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr rainfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr ashfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velXGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velYGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velZGrid = openvdb::FloatGrid::create();
            cloudGrid->setName("cloudDensityGrid");
            grids.push_back(cloudGrid);
            cloudTypeGrid->setName("cloudTypeGrid");
            grids.push_back(cloudTypeGrid);
            particulateGrid->setName("particulateDensityGrid");
            grids.push_back(particulateGrid);
            rainfallGrid->setName("rainfallDensityGrid");
            grids.push_back(rainfallGrid);
            ashfallGrid->setName("ashfallDensityGrid");
            grids.push_back(ashfallGrid);
            velXGrid->setName("velXGrid");
            grids.push_back(velXGrid);
            velYGrid->setName("velYGrid");
            grids.push_back(velYGrid);
            velZGrid->setName("velZGrid");
            grids.push_back(velZGrid);

            auto cloud = cloudGrid->getAccessor();
            auto cType = cloudTypeGrid->getAccessor();
            auto parts = particulateGrid->getAccessor();
            auto rain = rainfallGrid->getAccessor();
            auto ash = ashfallGrid->getAccessor();
            auto velX = velXGrid->getAccessor();
            auto velY = velYGrid->getAccessor();
            auto velZ = velZGrid->getAccessor();
            int cellsPerHeight = wld->getAirLayer(0)->getObstacles().getY();
            int cellsPerWidth = wld->getAirLayer(0)->getObstacles().getX();
            auto cellWidth = wld->getCellSize();
            auto layers = wld->getAirLayers();
            auto convectLayers = wld->getConvectionLayers();
            double h, cH;

            for (int k = 0; k < layers.size(); ++k){
                h = layers[k]->getHeight();
                cloudGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                cloudGrid->insertMeta("cloudThreshold", openvdb::DoubleMetadata(cloudThres));
                cloudTypeGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                particulateGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velXGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velYGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                if (k < convectLayers.size()){
                    cH = convectLayers[k]->getLayerHeight();
                    rainfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    velZGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    ashfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                }
                for (int i = 0; i < cellsPerHeight; ++i){
                    for (int j = 0; j < cellsPerWidth; ++j){
//                        auto xy = wld->getWorldXYCoords(i, j);
                        openvdb::Coord ZYX(k, cellsPerHeight - 1 - i, j);
                        cloud.setValue(ZYX, layers[k]->getCondensedWater(i, j));
                        cType.setValue(ZYX, layers[k]->getClouds(i, j));
                        parts.setValue(ZYX, layers[k]->getParticulates(i, j));
                        velX.setValue(ZYX, layers[k]->getVelocityPhi(i, j));
                        velY.setValue(ZYX, -(layers[k]->getVelocityTheta(i, j)));
                        if (k < convectLayers.size()){
                            rain.setValue(ZYX, convectLayers[k]->getRainfall(i, j));
                            velZ.setValue(ZYX, convectLayers[k]->getVerticalVelocity(i, j));
                            ash.setValue(ZYX, convectLayers[k]->getAshfall(i, j));
                        }
                    }
                }
            }
            for (auto x : grids){
                x->insertMeta("secondsSinceStart", openvdb::DoubleMetadata(sim_time));
                x->insertMeta("cellsPerVerticalSide", openvdb::Int32Metadata(cellsPerHeight));
                x->insertMeta("cellsPerHorizontalSide", openvdb::Int32Metadata(cellsPerWidth));
                x->insertMeta("cellSize", openvdb::DoubleMetadata(cellWidth));
            }

            std::stringstream fG;
            fG << outputDirectory << "/Weather_At_Frame_" << step << ".vdb";
            try {
                std::cerr << "Writing weather data to " << fG.str() << "." << std::endl;
                openvdb::io::File f(fG.str());
//                std::cerr << "File " << fG.str() << " has compression " << f.compression() << std::endl;
//                f.setCompression(openvdb::io::COMPRESS_ZIP | openvdb::io::COMPRESS_BLOSC);
                f.write(grids);
                f.close();
//                double i = 0;
//                for (auto x : grids){
//                    i += x->memUsage();
//                }

                for (auto x : grids)
                    x->clear();

//                double j = 0;
//                for (auto x : grids){
//                    j += x->memUsage();
//                }
//                std::cerr << "Mem usage before clear: " << i << std::endl;
//                std::cerr << "Mem usage after clear: " << j << std::endl;
                return 0;
            } catch (std::exception& e){
                std::cerr << "\033[1;31mError! Can't write weather data to file " << fG.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mWeather data not written.\033[0m" << std::endl;
                for (auto x : grids){
                    x->clear();
                }
                return -1;
            }
        }

        template<typename T, typename TT, typename V, typename VV>
        int vdbExporter<T, TT, V, VV>::write(const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& airLayers, const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& convectionLayers, std::string outputDirectory, int step, double sim_time, V cloudThres){
            openvdb::GridPtrVec grids;
            openvdb::FloatGrid::Ptr cloudGrid = openvdb::FloatGrid::create();
            openvdb::Int32Grid::Ptr cloudTypeGrid = openvdb::Int32Grid::create();
            openvdb::FloatGrid::Ptr particulateGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr rainfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr ashfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velXGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velYGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velZGrid = openvdb::FloatGrid::create();
            cloudGrid->setName("cloudDensityGrid");
            grids.push_back(cloudGrid);
            cloudTypeGrid->setName("cloudTypeGrid");
            grids.push_back(cloudTypeGrid);
            particulateGrid->setName("particulateDensityGrid");
            grids.push_back(particulateGrid);
            rainfallGrid->setName("rainfallDensityGrid");
            grids.push_back(rainfallGrid);
            ashfallGrid->setName("ashfallDensityGrid");
            grids.push_back(ashfallGrid);
            velXGrid->setName("velXGrid");
            grids.push_back(velXGrid);
            velYGrid->setName("velYGrid");
            grids.push_back(velYGrid);
            velZGrid->setName("velZGrid");
            grids.push_back(velZGrid);

            auto cloud = cloudGrid->getAccessor();
            auto cType = cloudTypeGrid->getAccessor();
            auto parts = particulateGrid->getAccessor();
            auto rain = rainfallGrid->getAccessor();
            auto ash = ashfallGrid->getAccessor();
            auto velX = velXGrid->getAccessor();
            auto velY = velYGrid->getAccessor();
            auto velZ = velZGrid->getAccessor();
            int cellsPerSide = airLayers.at(0)->getObstacles().getX();
            auto cellWidth = airLayers.at(0)->getObstacles().cellSizeX();
            V h, cH;
            for (int k = 0; k < airLayers.size(); ++k){
                h = airLayers[k]->getHeight();
                cloudGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                cloudGrid->insertMeta("cloudThreshold", openvdb::DoubleMetadata(cloudThres));
                cloudTypeGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                particulateGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velXGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velYGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                if (k < convectionLayers.size()){
                    cH = convectionLayers[k]->getLayerHeight();
                    velZGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    rainfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    ashfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                }
                for (int i = 0; i < cellsPerSide; ++i){
                    for (int j = 0; j < cellsPerSide; ++j){
                        //                        auto xy = wld->getWorldXYCoords(i, j);
                        openvdb::Coord ZYX(k, cellsPerSide - 1 - i, j);
                        cloud.setValue(ZYX, airLayers[k]->getCondensedWater(i, j));
                        cType.setValue(ZYX, airLayers[k]->getClouds(i, j));
                        parts.setValue(ZYX, airLayers[k]->getParticulates(i, j));
                        velX.setValue(ZYX, airLayers[k]->getVelocityPhi(i, j));
                        velY.setValue(ZYX, -(airLayers[k]->getVelocityTheta(i, j)));
                        if (k < convectionLayers.size()){
                            rain.setValue(ZYX, convectionLayers[k]->getRainfall(i, j));
                            velZ.setValue(ZYX, convectionLayers[k]->getVerticalVelocity(i, j));
                            ash.setValue(ZYX, convectionLayers[k]->getAshfall(i, j));
                        }
                    }
                }
            }
            for (auto x : grids){
                x->insertMeta("secondsSinceStart", openvdb::DoubleMetadata(sim_time));
                x->insertMeta("cellsPerSide", openvdb::Int32Metadata(cellsPerSide));
                x->insertMeta("cellSize", openvdb::DoubleMetadata(cellWidth));
            }


            std::stringstream fG;
            fG << outputDirectory << "/Weather_At_Frame_" << step << ".vdb";
            try {
                std::cerr << "Writing weather data to " << fG.str() << "." << std::endl;
                openvdb::io::File f(fG.str());
                f.write(grids);
                f.close();

                for (auto x : grids)
                    x->clear();
                return 0;
            } catch (std::exception& e){
                std::cerr << "\033[1;31mError! Can't write weather data to file " << fG.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mWeather data not written.\033[0m" << std::endl;
                for (auto x : grids){
                    x->clear();
                }
                return -1;
            }
        }

        template<>
        int vdbExporter<PWM::PWMDataStructure::flatStaggeredGrid<double>, PWM::PWMDataStructure::flatStaggeredGrid<std::string>, double, std::string>::write(const std::vector<std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>>& airLayers, const std::vector<std::shared_ptr<PWM::Model::convectionLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>>& convectionLayers, std::string outputDirectory, int step, double sim_time, double cloudThres){
            openvdb::GridPtrVec grids;
            openvdb::FloatGrid::Ptr cloudGrid = openvdb::FloatGrid::create();
            openvdb::Int32Grid::Ptr cloudTypeGrid = openvdb::Int32Grid::create();
            openvdb::FloatGrid::Ptr particulateGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr rainfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr ashfallGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velXGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velYGrid = openvdb::FloatGrid::create();
            openvdb::FloatGrid::Ptr velZGrid = openvdb::FloatGrid::create();
            cloudGrid->setName("cloudDensityGrid");
            grids.push_back(cloudGrid);
            cloudTypeGrid->setName("cloudTypeGrid");
            grids.push_back(cloudTypeGrid);
            particulateGrid->setName("particulateDensityGrid");
            grids.push_back(particulateGrid);
            rainfallGrid->setName("rainfallDensityGrid");
            grids.push_back(rainfallGrid);
            ashfallGrid->setName("ashfallDensityGrid");
            grids.push_back(ashfallGrid);
            velXGrid->setName("velXGrid");
            grids.push_back(velXGrid);
            velYGrid->setName("velYGrid");
            grids.push_back(velYGrid);
            velZGrid->setName("velZGrid");
            grids.push_back(velZGrid);

            auto cloud = cloudGrid->getAccessor();
            auto cType = cloudTypeGrid->getAccessor();
            auto parts = particulateGrid->getAccessor();
            auto rain = rainfallGrid->getAccessor();
            auto ash = ashfallGrid->getAccessor();
            auto velX = velXGrid->getAccessor();
            auto velY = velYGrid->getAccessor();
            auto velZ = velZGrid->getAccessor();
            int cellsPerSide = airLayers.at(0)->getObstacles().getX();
            auto cellWidth = airLayers.at(0)->getObstacles().cellSizeX();
            double h, cH;
            for (int k = 0; k < airLayers.size(); ++k){
                h = airLayers[k]->getHeight();
                cloudGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                cloudGrid->insertMeta("cloudThreshold", openvdb::DoubleMetadata(cloudThres));
                cloudTypeGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                particulateGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velXGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                velYGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(h));
                if (k < convectionLayers.size()){
                    cH = convectionLayers[k]->getLayerHeight();
                    velZGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    rainfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                    ashfallGrid->insertMeta("heightAt_" + std::to_string(k), openvdb::DoubleMetadata(cH));
                }
                for (int i = 0; i < cellsPerSide; ++i){
                    for (int j = 0; j < cellsPerSide; ++j){
                        //                        auto xy = wld->getWorldXYCoords(i, j);
                        openvdb::Coord ZYX(k, cellsPerSide - 1 - i, j);
                        cloud.setValue(ZYX, airLayers[k]->getCondensedWater(i, j));
                        cType.setValue(ZYX, airLayers[k]->getClouds(i, j));
                        parts.setValue(ZYX, airLayers[k]->getParticulates(i, j));
                        velX.setValue(ZYX, airLayers[k]->getVelocityPhi(i, j));
                        velY.setValue(ZYX, -(airLayers[k]->getVelocityTheta(i, j)));
                        if (k < convectionLayers.size()){
                            rain.setValue(ZYX, convectionLayers[k]->getRainfall(i, j));
                            velZ.setValue(ZYX, convectionLayers[k]->getVerticalVelocity(i, j));
                            ash.setValue(ZYX, convectionLayers[k]->getAshfall(i, j));
                        }
                    }
                }
            }
            for (auto x : grids){
                x->insertMeta("secondsSinceStart", openvdb::DoubleMetadata(sim_time));
                x->insertMeta("cellsPerSide", openvdb::Int32Metadata(cellsPerSide));
                x->insertMeta("cellSize", openvdb::DoubleMetadata(cellWidth));
            }


            std::stringstream fG;
            fG << outputDirectory << "/Weather_At_Frame_" << step << ".vdb";
            try {
                std::cerr << "Writing weather data to " << fG.str() << "." << std::endl;
                openvdb::io::File f(fG.str());
                f.write(grids);
                f.close();

                for (auto x : grids)
                    x->clear();
                return 0;
            } catch (std::exception& e){
                std::cerr << "\033[1;31mError! Can't write weather data to file " << fG.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mWeather data not written.\033[0m" << std::endl;
                for (auto x : grids){
                    x->clear();
                }
                return -1;
            }
        }
    }
}

#endif // PWM_UTILS_VDBEXPORTER_H
