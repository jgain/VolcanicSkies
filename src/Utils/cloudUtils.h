#ifndef PWM_UTILS_CLOUDS
#define PWM_UTILS_CLOUDS

#include <array>
#include <chrono>
#include "cloudType.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "vizUtils.h"
#include "world.h"

#define cumulusThreshold 0.5
#define stratusThreshold 0.15

#define nimboRange 2000
#define cumulonimbusRange 1500

#define noPrefixHeight 2000
#define altoHeight 6000

#define tropopause 10000
#define waterLimit 10

#define kMax 15
#define kUpliftTol 0.1
#define kWaterTol 0.1
#define kAngleTol 0.31415926535
#define kDiffTol 0.00001
#define kMaxIter 50
#define kTestTol 0.25
#define kReps 5
#define kNumSeeds 30

namespace PWM{
    namespace Utils{
        template <typename V>
        struct cloud{
            cloudType cType;
            V minAltitude;
            V maxAltitude;
        };

        template<typename V>
    inline cloudType cloudTaxonomy(V uplift, V base, V top){
            cloudType cType;
            cType = cloudType::EMPTY;

            if (uplift < stratusThreshold){
                if (top - base > nimboRange)
                    cType = cloudType::NIMBOSTRATUS;
                else{
                    if (top < noPrefixHeight)
                        cType = cloudType::STRATUS;
                    else if (top < altoHeight)
                        cType = cloudType::ALTOSTRATUS;
                    else
                        cType = cloudType::CIRROSTRATUS;
                }
            }
            else if (uplift > cumulusThreshold){
                if (top - base > cumulonimbusRange)
                    cType = cloudType::CUMULONIMBUS;
                else{
                    if (top < noPrefixHeight)
                        cType = cloudType::CUMULUS;
                    else if (top < altoHeight)
                        cType = cloudType::ALTOCUMULUS;
                    else
                        cType = cloudType::CIRROCUMULUS;
                }
            }
            else{
                if (top < noPrefixHeight)
                    cType = cloudType::STRATOCUMULUS;
                else if (top < altoHeight)
                    cType = cloudType::NIMBOSTRATUS;
                else
                    cType = cloudType::CIRRUS;
            }
            return cType;
        }

/*        template<typename T, typename V>
//        void writeLayerType(std::string outputFile, T* cloudMap, V base, V thickness, V radius){
//            std::ofstream outFile;
//            outFile.open((char *) outputFile.c_str(), std::ios_base::out);
//            if (outFile.is_open()){
//                outFile << base << " " << thickness << std::endl;
//                outFile << radius << std::endl;
//                outFile << cloudMap->size() << " cells" << std::endl;

//                for(int i = 0; i < cloudMap->size(); ++i)
//                    outFile << cloudMap->getData(i) << " ";
//            }
//            else
//                std::cerr << "Error in writing cloud type (cloudUtils::writeLayerType()) -- can't open " << outputFile << std::endl;
//        }

//        template<typename T, typename TT, typename V, typename VV>
//        void writeVelocityField(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputFile, int layer){
//            std::ofstream outFile;
//            outFile.open((char *) outputFile.c_str(), std::ios_base::out);
//            if (outFile.is_open()){
//                outFile << wld->getAirLayer(layer)->getObstacles().size() << " cells" << std::endl;

//                for (int i = 0; i < wld->getAirLayer(layer)->getObstacles().size(); ++i)
//                    outFile << wld->getAirLayer(layer)->getVelocityTheta(i) << " " << wld->getAirLayer(layer)->getVelocityPhi(i) << " ";
//            }
//            else
//                std::cerr << "Error in writing velocity field (cloudUtils::writeVelocityField()) -- can't open " << outputFile << std::endl;
//        }

//        template<typename T, typename TT, typename V, typename VV>
//        void writeCloud(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputFile, T* densityMap, T* cloudMap, V base, V top, int cloudIndex, int layer){
//            T submap(*densityMap);
//            for(int i = 0; i < submap.size(); ++i){
//                V val = int(std::floor(255.0 * submap.getData(i)));
//                val = std::max(std::min(val, 255.0), 0.0);
//                submap.setData(i, val);
//            }
//            std::string dOutputFile = outputFile + "_d_" + std::to_string(cloudIndex) + ".png";

//            std::ofstream dOutput(dOutputFile);
//            dOutput << submap.print();
//            dOutput.close();

//            writeVelocityField<T, TT, V, VV>(wld, (outputFile + "_v_" + std::to_string(cloudIndex) + ".png"), layer);
//            writeLayerType<T, V>((outputFile + "_t_" + std::to_string(cloudIndex) + ".png"), cloudMap, base, top - base, wld->getPlanet()->getRadius());
//        }

//        template<typename T, typename TT, typename V, typename VV>
//        void segmentCloud(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputFile, T* mask, int cloudIndex, int layer, V base, V top){
//            V max = 0.0f;
//            T densityMap(wld->getAirLayer(layer)->getObstacles());
//            T cloudMap(wld->getAirLayer(layer)->getObstacles());
//            densityMap = 0;
//            cloudMap = 0;
//            V uplift;
//            cloudType cType;

//            for (int i = 0; i < wld->getAirLayer(layer)->getObstacles().size(); ++i){
//                if (mask->getData(i) == cloudIndex){
//                    if (layer > 0)
//                        uplift = std::abs(wld->getConvectionLayer(layer - 1)->getVerticalVelocity(i));
//                    else
//                        uplift = 0.0f;

//                    cType = cloudTaxonomy<V>(uplift, base, top);
//                    cloudMap.setData(i, (int) cType);
//                }
//                else
//                    cloudMap.setData(i, (int) cloudType::EMPTY);
                
//                if (mask->getData(i) == cloudIndex){
//                    densityMap.setData(i, wld->getAirLayer(layer)->getCondensedWater(i));
//                    if (densityMap.getData(i) > max)
//                        max = densityMap.getData(i);
//                }
//                else
//                    densityMap.setData(i, 0);
//            }
//            #pragma omp parallel for
//                for (int i = 0; i < wld->getAirLayer(layer)->getObstacles().size(); ++i){
//                    densityMap.setData(i, std::min(1.0, densityMap.getData(i) / waterLimit));
//                }
//            #pragma omp barrier

//            writeCloud(wld, outputFile, &densityMap, &cloudMap, base, top, cloudIndex, layer);
//        }*/

        template<typename V>
    inline V kMeansDist(std::array<V, 4> f1, std::array<V, 4> f2){
            V diff = 0.0f;
            for(int e = 0; e < 4; ++e)
                diff += std::pow(f1[e] - f2[e], 2);
            
            return diff;
        }

        template<typename V>
    inline bool kMeansClose(int k, std::vector<std::array<V, 4>>& clusters, V tol){
            for (int i = 0; i < k; ++i)
                for (int j = i + 1; j < k; ++j)
                    if (std::sqrt(kMeansDist(clusters[i], clusters[j])) < tol)
                        return true;
            return false;
        }

        template<typename V>
    inline V kMeans(std::vector<std::array<V, 4>>& features, int k, std::vector<int>& assigns, std::vector<std::array<V, 4>>& clusters){
            std::vector<int> kIndices;
            std::default_random_engine generator(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
            std::uniform_int_distribution distribution(std::uniform_int_distribution<int>(0, (int) features.size() - 1));
            std::vector<int> prevAssigns;
            V diff = 0.0f, bestDist;
            bool converged;
            std::vector<int> rndSeeds;

            clusters.clear();
            rndSeeds.assign(kNumSeeds, -1);

            if ((int) features.size() <= kNumSeeds){
                std::cerr << "Error in kMeans! Feature vector less than number of seeds in kMeans clustering. Returning default value." << std::endl;
                return diff;
            }

            for (int c = 0; c < k; ++c){
                int rnd, bIndex;
                for (int r = 0; r < kNumSeeds; ++r){
                    bool valid = false;
                    while(!valid){
                        rnd = distribution.operator()(generator);
                        valid = true;
                        for (int e = 0; e < c; ++e)
                            if (kIndices[e] == rnd){
                                valid = false;
                                break;
                            }
                    }
                    rndSeeds[r] = rnd;
                }

                bestDist = 0.0f;
                bIndex = rndSeeds[0];

                for (int r = 0; r < kNumSeeds; ++r){
                    V shortDist = std::numeric_limits<V>::max();
                    for (int i = 0; i < c; ++i){
                        V dist = std::sqrt(kMeansDist<V>(features[rndSeeds[r]], clusters[i]));
                        if (dist < shortDist)
                            shortDist = dist;
                    }

                    if (shortDist > bestDist){
                        bIndex = rndSeeds[r];
                        bestDist = shortDist;
                    }
                }

                kIndices.push_back(bIndex);
                clusters.push_back(features[bIndex]);
            }

            assigns.clear();
            assigns.assign(features.size(), 0);
            prevAssigns.assign(features.size(), 0);

            converged = false;
            int j = 0;
            while (!converged && j < kMaxIter){
                for (int f = 0; f < (int) features.size(); ++f){
                    V bestDist = std::numeric_limits<V>::max();
                    V currDist;
                    int cIndex = 0;

                    for (int c = 0; c < k; ++c){
                        currDist = kMeansDist<V>(features[f], clusters[c]);
                        if (currDist < bestDist){
                            bestDist = currDist;
                            cIndex = c;
                        }
                    }

                    assigns[f] = cIndex;
                }

                for (int c = 0; c < k; ++c){
                    int kCount = 0;
                    for (int e = 0; e < 4; ++e)
                        clusters[e][c] = 0.0f;
                    
                    for (int f = 0; f < (int) features.size(); ++f){
                        if (assigns[f] == c){
                            for (int e = 0; e < 4; ++e)
                                clusters[c][e] += features[f][e];
                            ++kCount;
                        }
                    }
                    for (int e = 0; e < 4; ++e)
                        clusters[c][e] /= (V) kCount;
                }
                converged = true;
                for (int f = 0; f < (int) features.size(); ++f){
                    if (prevAssigns[f] != assigns[f])
                        converged = false;
                    prevAssigns[f] = assigns[f];
                }

                ++j;
            }

            for (int f = 0; f < features.size(); ++f)
                diff += kMeansDist<V>(features[f], clusters[assigns[f]]);
            
            diff /= (V) features.size();
            return diff;
        }

        template<typename V>
    inline void kMeansSearch(std::vector<std::array<V, 4>>& features, int& k, std::vector<int>& assigns, std::vector<std::array<V, 4>>& clusters){
            bool fin = false;
            V currDiff, bestDiff, globalDiff;
            std::vector<int> bestAssigns, currAssigns;
            std::vector<std::array<V, 4>> bestClusters, currClusters;

            globalDiff = std::numeric_limits<V>::max();

            for (int c = 1; c < kMax; ++c){
                bestDiff = std::numeric_limits<V>::max();
                for (int r = 0; r < kReps; ++r){
                    currDiff = kMeans<V>(features, c, currAssigns, currClusters);
                    if (currDiff < bestDiff && !kMeansClose<V>(c, currClusters, 0.06f)){
                        bestDiff = currDiff;
                        bestAssigns = currAssigns;
                        bestClusters = currClusters;
                    }
                }

                if (bestDiff < globalDiff){
                    k = c;
                    assigns = bestAssigns;
                    clusters = bestClusters;
                }
            }
        }

        template <typename T, typename TT, typename V, typename VV>
    inline void basicCloudExtract(const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputDirectory, int step, V minCondensationThreshold){
//            std::cerr << "Starting cloud extraction" << std::endl;
//            std::vector<std::shared_ptr<T>> masks;
            size_t cellCount = wld->getAirLayer(0)->getObstacles().size();
            size_t layerCount = wld->getAirLayers().size();

            V noCloudVal = 0;
            V cloudUnassignedVal = -1;
            
            for (auto l : wld->getAirLayers()){
//                auto lMask = std::make_shared<T>(l->getObstacles());

                 #pragma omp parallel for
                    for (int i = 0; i < cellCount; ++i){
                        if (l->getObstacles(i)){
                            l->setClouds(i, noCloudVal);
                            continue;
                        }
                        V condensed = l->getCondensedWater(i);
                        if (condensed <= minCondensationThreshold){
                            l->setClouds(i, noCloudVal);
                        }
                        else{
                            l->setClouds(i, cloudUnassignedVal);
                        }
                    }
                 #pragma omp barrier
//                masks.push_back(lMask);
//                std::cout << "Average temperature for layer " << masks.size() << ": " << l->getMeanTemperature() << " K." << std::endl;
            }
            
            for (int i = 0; i < layerCount; ++i){
                #pragma omp parallel for
                for (int j = 0; j < cellCount; ++j){
                    if (wld->getAirLayer(i)->getClouds(j) != noCloudVal){
                        V uplift = 0;
                        if (i == 0)
                            uplift = wld->getConvectionLayer(i)->getVerticalVelocity(j);
                        else if (i == layerCount - 1)
                            uplift = wld->getConvectionLayer(i - 1)->getVerticalVelocity(j);
                        else
                            uplift = 0.5 * (wld->getConvectionLayer(i - 1)->getVerticalVelocity(j) + wld->getConvectionLayer(i)->getVerticalVelocity(j));
                        V base = wld->getAirLayer(i)->getLayerBot();
                        V top = wld->getAirLayer(i)->getLayerTop();
                        for (int k = i + 1; k < layerCount; ++k){
                            if (wld->getAirLayer(k)->getClouds(j) != noCloudVal){
                                top = wld->getAirLayer(k)->getLayerTop();
                            }
                            else
                                break;
                        }
                        for (int k = i - 1; k >= 0; --k){
                            if (wld->getAirLayer(k)->getClouds(j) != noCloudVal){
                                base = wld->getAirLayer(k)->getLayerBot();
                            }
                            else
                                break;
                        }
                        auto cType = cloudTaxonomy(uplift, base, top);
                        wld->getAirLayer(i)->setClouds(j, static_cast<V>(cType));
                    }
                }
//                #pragma omp barrier
//                writeVelocityField<T, TT, V, VV>(wld, (outputFile + "_v_" + std::to_string(cloudIndex) + ".png"), i - 1);
//                writeLayerType((outputFile + "_t_.png"), masks[i - 1], wld->getAirLayer(i)->getHeight(), 0.0, wld->getPlanet()->getRadius());
            }


//            for (int i = 0; i < layerCount; ++i){
//                writeCloudImages<V>(outputDirectory, step, i, wld->getAirLayer(i)->getHeight(), wld->getAirLayer(i)->getClouds());
//                wld->getAirLayer(i)->getClouds().copy(*(masks[i]));
//            }
        }

        template <typename T, typename V>
        inline void basicCloudExtract(const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& airLayers, const std::vector<std::shared_ptr<PWM::Model::convectionLayer<T, V>>>& convectLayers, V minCondensationThreshold){
            size_t cellCount = airLayers.at(0)->getObstacles().size();
            size_t layerCount = airLayers.size();

            V noCloudVal = 0;
            V cloudUnassignedVal = -1;

            for (auto l : airLayers){
                //                auto lMask = std::make_shared<T>(l->getObstacles());

                #pragma omp parallel for
                for (int i = 0; i < cellCount; ++i){
                    if (l->getObstacles(i)){
                        l->setClouds(i, noCloudVal);
                        continue;
                    }
                    V condensed = l->getCondensedWater(i);
                    if (condensed <= minCondensationThreshold){
                        l->setClouds(i, noCloudVal);
                    }
                    else{
                        l->setClouds(i, cloudUnassignedVal);
                    }
                }
                #pragma omp barrier
            }

            for (int i = 0; i < layerCount; ++i){
                #pragma omp parallel for
                for (int j = 0; j < cellCount; ++j){
                    if (airLayers.at(i)->getClouds(j) != noCloudVal){
                        V uplift = 0;
                        if (i == 0)
                            uplift = convectLayers.at(i)->getVerticalVelocity(j);
                        else if (i == layerCount - 1)
                            uplift = convectLayers.at(i - 1)->getVerticalVelocity(j);
                        else
                            uplift = 0.5 * (convectLayers.at(i - 1)->getVerticalVelocity(j) + convectLayers.at(i)->getVerticalVelocity(j));
                        V base = airLayers.at(i)->getLayerBot();
                        V top = airLayers.at(i)->getLayerTop();
                        for (int k = i + 1; k < layerCount; ++k){
                            if (airLayers.at(k)->getClouds(j) != noCloudVal){
                                top = airLayers.at(k)->getLayerTop();
                            }
                            else
                                break;
                        }
                        for (int k = i - 1; k >= 0; --k){
                            if (airLayers.at(k)->getClouds(j) != noCloudVal){
                                base = airLayers.at(k)->getLayerBot();
                            }
                            else
                                break;
                        }
                        auto cType = cloudTaxonomy(uplift, base, top);
                        airLayers.at(i)->setClouds(j, static_cast<V>(cType));
                    }
                }
                #pragma omp barrier
            }
        }

//        template <typename T, typename TT, typename V, typename VV>
//        void extractClouds(std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld, std::string outputFile, V minCondensationThreshold){
//            // Define data structures - per layer, a mask of cloud indices
//            // Per cloud index, type, density, min elev, max elev, contained in struct
//            std::cerr << "Starting cloud extraction" << std::endl;
//            std::vector<T *> masks;
//            std::vector<cloud<V>> clouds;

//            int cloudIndex = 1;
//            V noCloudVal = -1;
            
//            for (auto l : wld->getAirLayers()){
//                T* lMask = new T;
//                *lMask = l->getObstacles();

//                // #pragma omp parallel for
//                    for (int i = 0; i < lMask->size(); ++i){
//                        V condensed = l->getCondensedWater(i);
//                        if (condensed <= minCondensationThreshold)
//                            lMask->setData(i, noCloudVal);
//                    }
//                // #pragma omp barrier

//                masks.push_back(lMask);
//            }

//            std::cerr << "Starting cloud locating" << std::endl;

//            for (int i = 1; i < wld->getAirLayers().size(); ++i){
//                std::vector<std::array<V, 4>> features, clusters;
//                std::vector<int> assigns;
//                bool cumulusFound;
//                V lowWater, medWater, highWater, uplift, lowHeight, medHeight, highHeight, base, top, w;
//                int k, count;
//                cloudType cType;

//                features.clear();
//                clusters.clear();
//                assigns.clear();

//                std::cerr << "Starting cumulus found section for layer " << i << std::endl;
//                cumulusFound = false;
//                for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                    if (std::abs(wld->getConvectionLayer(i - 1)->getVerticalVelocity(j) > cumulusFound)){
//                        masks[i - 1]->setData(j, cloudIndex);
//                        cumulusFound = true;
//                        break;
//                    }
//                }
//                // #pragma omp parallel for
//                    for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                        if (masks[i - 1]->getData(j) != noCloudVal)
//                            masks[i - 1]->setData(j, cloudIndex);
//                    }
//                // #pragma omp barrier
//                //cumulusFound = true;
//                std::cerr << "performing cumulus found section for layer " << i << std::endl;
//                if(cumulusFound){
//                    lowWater = 0.0f;
//                    medWater = 0.0f;
//                    highWater = 0.0f;
//                    uplift = 0.0f;

//                    // Can I use a conditional reduction here?
//                    // #pragma omp parallel for reduction(+: lowWater, medWater, highWater, uplift, count)
//                        for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                            if (masks[i - 1]->getData(j) == noCloudVal){
//                                lowWater += wld->getAirLayer(i - 1)->getCondensedWater(j);
//                                medWater += wld->getAirLayer(i)->getCondensedWater(j);
//                                if (i + 1 < wld->getAirLayers().size())
//                                    highWater += wld->getAirLayer(i + 1)->getCondensedWater(j);
//                                uplift += wld->getConvectionLayer(i - 1)->getVerticalVelocity(j);
//                                count += 1;
//                            }
//                        }
//                    // #pragma omp barrier

//                    lowWater /= (V) count;
//                    medWater /= (V) count;
//                    highWater /= (V) count;
//                    uplift /= (V) count;

//                    if (i > 0)
//                        lowHeight = wld->getAirLayer(i - 1)->getHeight();
//                    else
//                        lowHeight = 0.0;
                    
//                    medHeight = wld->getAirLayer(i)->getHeight();
//                    if (i + 1 < wld->getAirLayers().size())
//                        highHeight = wld->getAirLayer(i + 1)->getHeight();
//                    else
//                        highHeight = tropopause;
                    
//                    base = (lowHeight + medHeight) / 2.0;
//                    top = (medHeight + highHeight) / 2.0;
//                    segmentCloud(wld, outputFile, masks[i - 1], cloudIndex, i, base, top);
//                    ++cloudIndex;
//                }

//                std::cerr << "Starting feature assignment section for layer " << i << std::endl;
//                for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                    if (masks[i - 1]->getData(j) == 0){
//                        std::array<V, 4> f;
//                        f[0] = wld->getAirLayer(i)->getCondensedWater(j);
//                        f[1] = wld->getConvectionLayer(i - 1)->getVerticalVelocity(j);
//                        f[2] = wld->getAirLayer(i)->getVelocityTheta(i);
//                        f[3] = wld->getAirLayer(i)->getVelocityPhi(i);
//                        features.push_back(f);
//                    }
//                }
                
//                std::cerr << "Starting kMeansSearch section for layer " << i << std::endl;
//                kMeansSearch<V>(features, k, assigns, clusters);

//                for (int c = 0; c < clusters.size(); c++){
//                    int baseLayer, topLayer, currLayer;
//                    int f = 0;

//                    lowWater = 0.0f;
//                    count = 0;
//                    // Can I use a conditional reduction here? Doesn't look like it, f is a problem.
//                    // #pragma omp parallel for reduction(+: lowWater, count)
//                        for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                            if (masks[i - 1]->getData(j) == 0){
//                                if (assigns[f] == c){
//                                    lowWater += wld->getAirLayer(i - 1)->getCondensedWater(j);
//                                    ++count;
//                                }
//                                ++f;
//                            }
//                        }
//                    // #pragma omp barrier

//                    lowWater /= (V) count;
//                    highWater = clusters[c][0];
//                    lowHeight = wld->getAirLayer(i - 1)->getHeight();
//                    highHeight = wld->getAirLayer(i)->getHeight();

//                    if (lowWater > minCondensationThreshold){
//                        std::cerr << "Error in extractClouds(): missed cloud in layer " << i - 1 << std::endl;
//                        base = lowHeight;
//                    }
//                    else{
//                        w = (minCondensationThreshold - lowWater) / (highWater - lowWater);
//                        base = lowHeight + w * (highHeight - lowHeight);
//                    }

//                    currLayer = baseLayer;
//                    bool fin = false;
//                    while (!fin){
//                        lowWater = highWater;
//                        lowHeight = highHeight;
//                        highWater = 0.0;
//                        if (currLayer == wld->getAirLayers().size() - 1){
//                            highHeight = tropopause;
//                            fin = true;
//                        }
//                        else{
//                            highHeight = wld->getAirLayer(currLayer)->getHeight();
//                            f = 0;
//                            count = 0;
//                            // Can I use a conditional reduction here? Doesn't look like it, f is a problem.
//                            // #pragma omp parallel for reduction(+: lowWater, count)
//                                for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                                    if (masks[i - 1]->getData(j) == 0){
//                                        if (assigns[f] == c){
//                                            highWater += wld->getAirLayer(currLayer)->getCondensedWater(j);
//                                            ++count;
//                                        }
//                                        ++f;
//                                    }
//                                }
//                            // #pragma omp barrier
//                            highWater /= (V) count;
//                            if (highWater < minCondensationThreshold)
//                                fin = true;
//                        }
//                        ++currLayer;
//                    }
//                    topLayer = currLayer;

//                    w = (lowWater - minCondensationThreshold) / (lowWater - highWater);
//                    top = lowHeight + w * (highHeight - lowHeight);

//                    f = 0;
//                    // Can I use a conditional reduction here? Doesn't look like it, f is a problem.
//                    // #pragma omp parallel for reduction(+: lowWater, count)
//                        for (int j = 0; j < wld->getAirLayer(i)->getObstacles().size(); ++j){
//                            if (masks[i - 1]->getData(j) == 0 || masks[i - 1]->getData(j) < cloudIndex){
//                                if (assigns[f] == c){
//                                    for (int x = baseLayer; x < topLayer; ++x)
//                                        masks[x - 1]->setData(i, cloudIndex);
//                                }
//                                ++f;
//                            }
//                        }
//                    // #pragma omp barrier

//                    uplift = clusters[c][0];
//                    cType = cloudTaxonomy<V>(uplift, base, top);

//                    ++cloudIndex;
//                }
//            }

//            for (int i = 0; i < masks.size(); ++i)
//                delete masks[i];
//            std::cerr << "Cloud extraction cleanup complete." << std::endl;
//        }
    }
}

#endif //PWM_UTILS_CLOUDS
