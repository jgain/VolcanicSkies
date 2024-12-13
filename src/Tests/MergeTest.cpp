#include "cloudUtils.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include "PerlinNoise.hpp"
#include "pressureWave.h"
#include "smoke.hpp"
#include "skirt.hpp"
#include "terrain_structure.hpp"
#include "vdbExporter.h"
#include "world.h"
#include "worldEngine.h"
#include "timer.hpp"

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;
using namespace PWM::Engine;

int plumeToAtmoTransfers = 0;
std::pair<double, double> plumeLocation = std::make_pair(75000, 75000);
std::pair<double, double> BigPlumeLocation = std::make_pair(29500, 17500);
std::pair<double, double> smallPlumeLocation = std::make_pair(15000, 3000);

terrain_structure constructTestTerrStruct(std::string file, float cellWidth = 400.0f, float sceneWidth = 15000.0f){
    terrain_structure res(file, cellWidth, sceneWidth);
    return res;
}

/**
 * @brief constructBigTerrain
 * Constructs a full resolution terrain model from the given file.
 * @param file
 * @return
 */
terrain_structure constructBigTerrain(std::string file){
    return terrain_structure(file);
}

/**
 * @brief constructHalfTerrain
 * Constructs a half resolution terrain model from the given file.
 * @param file
 * @return
 */
terrain_structure constructHalfTerrain(std::string file){
    auto initTer = terrain_structure(file);

    auto finTer = terrain_structure();
    finTer.cell_size = initTer.cell_size * 2;
    finTer.field_size = initTer.field_size / 2;
    finTer.latitude = initTer.latitude;
    finTer.max_xyz = initTer.max_xyz;
    finTer.min_xyz = -initTer.max_xyz;

    for (int i = 0; i < initTer.field_size; i+=2){
        auto vec = std::vector<float>();
        auto vecN = std::vector<Eigen::Vector3f>();
        for (int j = 0; j < initTer.field_size; j+=2){
            vec.push_back(initTer.height_field[i][j]);
            vecN.push_back(initTer.normal_field[i][j]);
        }
        finTer.height_field.push_back(vec);
        finTer.normal_field.push_back(vecN);
    }
    return finTer;
}

terrain_structure constructFlatBigTerrain(float width, int cellCount){
    auto finTer = terrain_structure();
    finTer.cell_size = width * 2 / cellCount;
    finTer.field_size = cellCount;
    finTer.max_xyz = width;
    finTer.min_xyz = -width;
    finTer.latitude = 0;
    for (int i = 0; i < cellCount; ++i){
        auto vec = std::vector<float>();
        auto vecN = std::vector<Eigen::Vector3f>();
        for (int j = 0; j < cellCount; ++j){
            vec.push_back(0.f);
            vecN.push_back(Eigen::Vector3f(0, 0, 1));
        }
        finTer.height_field.push_back(vec);
        finTer.normal_field.push_back(vecN);
    }
    return finTer;
}

/**
 * @brief constructPartTerrain
 * Construct a lower resolution terrain model from the given file, with scale determined by the second parameter.
 * @param file
 * @param factor
 * @return
 */
terrain_structure constructPartTerrain(std::string file, uint factor = 2){
    auto initTer = terrain_structure(file);

    auto finTer = terrain_structure();
    finTer.cell_size = initTer.cell_size * factor;
    finTer.field_size = initTer.field_size / factor;
    finTer.latitude = initTer.latitude;
    finTer.max_xyz = initTer.max_xyz;
    finTer.min_xyz = -initTer.max_xyz;

    for (int i = 0; i < initTer.field_size; i+=factor){
        auto vec = std::vector<float>();
        auto vecN = std::vector<Eigen::Vector3f>();
        for (int j = 0; j < initTer.field_size; j+=factor){
            vec.push_back(initTer.height_field[i][j]);
            vecN.push_back(initTer.normal_field[i][j]);
        }
        finTer.height_field.push_back(vec);
        finTer.normal_field.push_back(vecN);
    }
    return finTer;
}

terrain_structure constructSmallBigTerrain(std::string file, int startX, int startY, int width, int height){
    auto initTer = terrain_structure(file);

    auto finTer = terrain_structure();
    finTer.cell_size = initTer.cell_size;
    finTer.field_size = width;
    finTer.latitude = initTer.latitude;
    finTer.max_xyz = finTer.cell_size * (finTer.field_size / 2);
    finTer.min_xyz = -finTer.max_xyz;

    for (int i = startY; i < startY + height; ++i){
        auto vec = std::vector<float>();
        for (int j = startX; j < startX + width; ++j){
            vec.push_back(initTer.height_field[i][j]);
        }
        finTer.height_field.push_back(vec);
    }
    return finTer;
}

bool writeLayerStruct(std::string folder, const std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& w){
    std::stringstream f;
    auto layers = w->getAirLayers();
    f << layers.size() << std::endl;
    for (auto l : layers){
        f << l->getHeight() << std::endl;
    }
    std::cerr << "Writing layer struct file" << std::endl;

    std::string fil = folder + "/layerstruct.txt";
    std::ofstream file(fil.c_str()); 
    file << f.str();
    file.close();
    return true;
}

bool writeTerrainAshVals(std::string folder, const std::shared_ptr<PWM::Model::terrain<dsType, dsSType, valType, valType2>>& ter, int stepCount, std::pair<double, double> plumeLoc){
    std::stringstream f;
    f << folder << "/Terrain_Ash_Step_" << stepCount << ".elv";
    PWM::Utils::writeTerrAshFile(f.str(), ter->getAsh(), plumeLoc);
    return true;
//    int x = ter->getAsh().getX(), y = ter->getAsh().getY();
//    f << x << " " << y << std::endl;
//    for (int i = 0; i < y; ++i){
//        for (int j = 0; j < x; ++j){
//            f << ter->getAsh(i * y + j) << " ";
//        }
//        f << std::endl;
//    }

//    std::stringstream fil;
//    std::ofstream file(fil.str().c_str());
//    file << f.str();
//    file.close();
}

void convertFromCoordinatesToGridIdx(float x, float y, unsigned int& i, unsigned int& j, const std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& w)
{
    float cellSize = w->getCellSize();
    i = (w->getWorldCellWidth() / 2) + x/cellSize;
    j = (w->getWorldCellWidth() / 2) + y/cellSize;
}

std::pair<int, int> convertFromCoordinatesToGridIdx(float x, float y, const std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& w){
    float cellSize = w->getCellSize();
    int i = (w->getWorldCellWidth() / 2) + x/cellSize;
    int j = (w->getWorldCellWidth() / 2) + y/cellSize;
    return std::make_pair(i, j);
}

void convertFromGridIdxToCoordinates(float& x, float& y, unsigned int i, unsigned int j, const std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& w)
{
    float cellSize = w->getCellSize();
    x = (i - (w->getWorldCellWidth() / 2)) * cellSize;
    y = (j - (w->getWorldCellWidth() / 2)) * cellSize;
}

void transferDataFromAtmToPlume(scene_model& plume, PWM::Engine::worldEngine<dsType, dsSType, valType, valType2>& b)
{
//    std::cout << "begin data transfer" << std::endl;
    // take location of plume with last slice emitted
    /*Eigen::Vector3f position = plume.smoke_layers[plume.smoke_layers.size()-1].center;
    float r = plume.smoke_layers[plume.smoke_layers.size()-1].r;
    // take cells around this position in the grid, compute the means and fill the plume structures
    unsigned int min_i=0, max_i=1, min_j=0, max_j=1;
    convertFromCoordinatesToGridIdx((float)position.x()-r, (float)position.y()-r, min_i, min_j, b.getWorldModel());
    convertFromCoordinatesToGridIdx((float)position.x()+r, (float)position.y()+r, max_i, max_j, b.getWorldModel());
    unsigned int min_i2=0, max_i2=1, min_j2=0, max_j2=1;
    convertFromCoordinatesToGridIdx((float)position.x()-2.*r, (float)position.y()-2.*r, min_i2, min_j2, b.getWorldModel());
    convertFromCoordinatesToGridIdx((float)position.x()+2.*r, (float)position.y()+2.*r, max_i2, max_j2, b.getWorldModel());*/
    for (unsigned int k = 0; k<b.getWorldModel()->getAirLayers().size(); k++)
    {
        // take closest plume slice
        unsigned int min_idx = 0;
        float min_dist = abs(plume.smoke_layers[0].center.z() - b.getWorldModel()->getAirLayer(k)->getHeight());
        for (unsigned int i_slice = 0; i_slice<plume.smoke_layers.size(); i_slice++)
        {
            float dist = abs(plume.smoke_layers[i_slice].center.z() - b.getWorldModel()->getAirLayer(k)->getHeight());
            if (dist < min_dist)
            {
                min_dist = dist;
                min_idx = i_slice;
            }
        }
        Eigen::Vector3f position = plume.smoke_layers[min_idx].center;
        float r = plume.smoke_layers[min_idx].r;
        // take cells around this position in the grid, compute the means and fill the plume structures
        unsigned int min_i=0, max_i=1, min_j=0, max_j=1;
        convertFromCoordinatesToGridIdx((float)position.x()-2.*r, (float)position.y()-2.*r, min_i, min_j, b.getWorldModel());
        convertFromCoordinatesToGridIdx((float)position.x()+2.*r, (float)position.y()+2.*r, max_i, max_j, b.getWorldModel());
        unsigned int min_i2=0, max_i2=1, min_j2=0, max_j2=1;
        convertFromCoordinatesToGridIdx((float)position.x()-3.*r, (float)position.y()-3.*r, min_i2, min_j2, b.getWorldModel());
        convertFromCoordinatesToGridIdx((float)position.x()+3.*r, (float)position.y()+3.*r, max_i2, max_j2, b.getWorldModel());

        unsigned int nb_cells = 0;
        float sum_temp = 0.f;
        float sum_pres = 0.f;
        float sum_moist = 0.f;
        Eigen::Vector2f sum_wind = Eigen::Vector2f(0.f,0.f);
        for (unsigned int i = min_i2; i<=max_i2; i++)
        {
            for (unsigned int j = min_j2; j<=max_j2; j++)
            {
                if ((i<min_i || i>max_i) && (j<min_j || j>max_j))
                {
                    sum_temp += b.getWorldModel()->getAirLayer(k)->getTemperature(i,j);
                    sum_pres += b.getWorldModel()->getAirLayer(k)->getPressure(i,j);
                    sum_moist += b.getWorldModel()->getAirLayer(k)->getMoisture(i,j);
                    sum_wind += Eigen::Vector2f(b.getWorldModel()->getAirLayer(k)->getVelocityTheta(i,j), b.getWorldModel()->getAirLayer(k)->getVelocityPhi(i,j));
                    nb_cells++;
                }
            }
        }

        float mean_temp = sum_temp/nb_cells;
        float mean_pres = sum_pres/nb_cells;
        float mean_moist = sum_moist/nb_cells;
        Eigen::Vector2f mean_wind = sum_wind/nb_cells;

//        if (std::abs(mean_temp) >= 1500){
//            std::cerr << "Mean_temp out of range value! Sum_temp = " << sum_temp << " K, nb_cells = " << nb_cells << " cells, mean_temp = " << mean_temp << " K." << std::endl;
//            assert(std::abs(mean_temp) < 1500);
//        }
//        if (std::abs(mean_pres) >= 1000000000){
//            std::cerr << "Mean_pres out of range value! Sum_pres = " << sum_pres << " Pa, nb_cells = " << nb_cells << " cells, mean_pres = " << mean_pres << " Pa." << std::endl;
//            assert(std::abs(mean_pres) < 100000000);
//        }
//        if (std::abs(mean_wind.x()) >= 100){
//            std::cerr << "Mean_wind.x() out of range value! Sum_wind.x() = " << sum_wind.x() << " m.s^-1, nb_cells = " << nb_cells << " cells, mean_wind.x() = " << mean_wind.x() << " m.s^-1." << std::endl;
//            assert(std::abs(mean_wind.x()) < 100);
//        }
//        if (std::abs(mean_wind.y()) >= 100){
//            std::cerr << "Mean_wind.y() out of range value! Sum_wind.y() = " << sum_wind.y() << " m.s^-1, nb_cells = " << nb_cells << " cells, mean_wind.y() = " << mean_wind.y() << " m.s^-1." << std::endl;
//            assert(std::abs(mean_wind.y()) < 100);
//        }

        plume.temperatures[k] = mean_temp;
        plume.pressures[k] = mean_pres;
        plume.moistures[k] = mean_moist;
        plume.winds[k] = wind_structure(mean_wind.x(), mean_wind.y(), true);
//        std::cerr << "Number of cells found for transfer: " << nb_cells << ", from min_i = (" << min_i << ", " << min_i2 << ") and min_j = (" << min_j << ", " << min_j2 << ")\nto max_i = (" << max_i << ", " << max_i2 << ") and max_j = (" << max_j << ", " << max_j2 << ")." << std::endl;
    }
//    std::cout << "end data transfer" << std::endl;
}

void transferDataFromPlumeToAtm(scene_model& plume, PWM::Engine::worldEngine<dsType, dsSType, valType, valType2>& b)
{
    // for each atm layer, take closest plume slice (if there are any close enough)
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        int plume_slice_idx = -1;
        float min_dist = 10000.f;
        for (unsigned int j = 0; j<plume.smoke_layers.size(); j++)
        {
            float dist = std::abs(b.getWorldModel()->getAirLayer(i)->getHeight() - plume.smoke_layers[j].center.z());
            if (dist < min_dist)
            {
                min_dist = dist;
                plume_slice_idx = j;
            }
        }
        // if too far make as if there is no plume at this place
        if (min_dist > plume.smoke_layers[plume_slice_idx].r)
            plume_slice_idx = -1;

        // if there is a plume slice transfer data to atm grids where the slice is
        if (plume_slice_idx != -1)
        {
            Eigen::Vector3f position = plume.smoke_layers[plume_slice_idx].center;
            float r = plume.smoke_layers[plume_slice_idx].r;

            unsigned int min_i=0, max_i=1, min_j=0, max_j=1;
            convertFromCoordinatesToGridIdx((float)position.x()-r, (float)position.y()-r, min_i, min_j, b.getWorldModel());
            convertFromCoordinatesToGridIdx((float)position.x()+r, (float)position.y()+r, max_i, max_j, b.getWorldModel());

            for (unsigned int k = min_i; k<=max_i; k++)
            {
                for (unsigned int l = min_j; l<=max_j; l++)
                {
                    b.getWorldModel()->getAirLayer(i)->setTemperature(k, l, plume.smoke_layers[plume_slice_idx].temperature);
                }
            }
        }
    }
}

int findAirLayerIndex(const std::vector<std::shared_ptr<PWM::Model::airLayer<dsType, valType>>>& layers, valType height){
    int r = -1;
    for (int i = 0; i < layers.size(); ++i){
        if (height < (layers[i]->getHeight() + (layers[i]->getThickness() / 2))){
            r = i;
            return r;
        }
    }
    return r;
}

void stripEdgeAsh(const std::vector<std::shared_ptr<PWM::Model::airLayer<dsType, valType>>>& layers, const std::vector<std::shared_ptr<PWM::Model::convectionLayer<dsType, valType>>>& convectLayers){
    valType zero = 0;
    for (auto l : layers){
        int iMax = l->getParticulates().getX();
        int jMax = l->getParticulates().getY();
        for (int i = 0; i < iMax; ++i){
            l->setParticulates(i, 0, zero);
            l->setParticulates(i, jMax - 1, zero);
        }
        for (int j = 1; j < jMax; ++j){
            l->setParticulates(0, j, zero);
            l->setParticulates(iMax - 1, j, zero);
        }
    }
}

void transferAshFromPlumeToAtm(scene_model& plume, const std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& wld){
    //For each free sphere, add some ash to atmosphere cells in position
    int topLayerIndex = wld->getAirLayers().size() - 1;
    valType cellWidth = wld->getCellSize();
    valType minHeight = wld->getAirLayer(0)->getLayerBot();
    for (int k = 0; k < plume.free_spheres.size(); ++k){
        if (plume.free_spheres[k].center.z() < minHeight)
            continue;
        int radius = std::round(plume.free_spheres[k].r / cellWidth);
        valType sphereRadius = plume.free_spheres[k].r;
        valType fsAshLoss = plume.free_spheres[k].density_loss_last_dt;
        if (fsAshLoss <= 0 || sphereRadius <= 0)
            continue;
        auto layerLoc = plume.free_spheres[k].center;
        int layerIndex = findAirLayerIndex(wld->getAirLayers(), layerLoc.z());
        if (layerIndex == -1)//-1 means that the smoke layer has passed above the highest atmosphere layer
            layerIndex = topLayerIndex;
        // Easy enough to avoid by increasing the number of layers and atmotop, though the model behaviour becomes less accurate above 10 km

        //std::cerr << "Ash added to airlayer " << layerIndex << ", " << fsAshLoss << " kg per cubic metre, at time " << plume.t_sum << " seconds and height " << layerLoc.z() << " metres." << std::endl;
        auto layerCentre = convertFromCoordinatesToGridIdx((float) layerLoc.x(), (float) layerLoc.y(), wld);
        valType layerMassAdd = 0;
        for (int i = layerCentre.first - radius; i < layerCentre.first + radius; ++i){
            for (int j = layerCentre.second - radius; j < layerCentre.second + radius; ++j){ // JG Bug?
                valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
                if (dist <= radius){
                    auto newAsh = wld->getAirLayer(layerIndex)->getParticulates(i, j) + fsAshLoss;
                    wld->getAirLayer(layerIndex)->setParticulates(i, j, newAsh);
                }
            }
        }
        plume.free_spheres[k].density_loss_last_dt = 0;
    }
}

void transferDataFromPlumeToAtmAlt(scene_model& plume, const std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& wld){
    ++plumeToAtmoTransfers;

    int layerCount = wld->getAirLayers().size();

    // new: instead of considering each plume slice, take all atm layers and then the closest plume slice
    for (unsigned int k = 0; k<wld->getAirLayers().size(); k++)
    {
        // take closest plume slice
        unsigned int min_idx = 0;
        auto airHeight = wld->getAirLayer(k)->getHeight();
        float min_dist = abs(plume.smoke_layers[0].center.z() - airHeight);
        for (unsigned int i_slice = 0; i_slice<plume.smoke_layers.size(); i_slice++)
        {
            float dist = abs(plume.smoke_layers[i_slice].center.z() - wld->getAirLayer(k)->getHeight());
            if (dist < min_dist)
            {
                min_dist = dist;
                min_idx = i_slice;
            }
        }

        // find all spheres close enough to the layer
        std::vector<unsigned int> close_spheres;
        for (unsigned int j = 0; j<plume.free_spheres.size(); j++)
        {
            float height_diff = abs(plume.free_spheres[j].center.z() - wld->getAirLayer(k)->getHeight());
            if (height_diff < plume.free_spheres[j].r)
            {
                close_spheres.push_back(j);
            }
        }
        // if a cell is located in these spheres, it takes the temp/moist value of the closest plume slice

        int layerIndex = (int)k;
        int radius = std::round(plume.smoke_layers[min_idx].r / wld->getCellSize());
        if (plume.smoke_layers[min_idx].r < 0){
            std::cerr << "Radius is negative (" << plume.smoke_layers[min_idx].r << " m)!" << std::endl;
            continue;
        }
        else if (radius > wld->getWorldXSize()){
            std::cerr << "Radius is larger than the entire scene!" << std::endl;
            continue;
        }
        valType sLTemp = plume.smoke_layers[min_idx].temperature;
        valType sLMoist = plume.smoke_layers[min_idx].moisture;
        valType sLVertVel = plume.smoke_layers[min_idx].v.z();
        valType velModFactor = std::max(0.0, std::min(100.0, 140 - sLVertVel));
        auto layerLoc = plume.smoke_layers[min_idx].center;
        auto layerCentre = convertFromCoordinatesToGridIdx((float) layerLoc.x(), (float) layerLoc.y(), wld);

        for (int i = layerCentre.first - 2*radius; i < layerCentre.first + 2*radius; ++i){
            for (int j = layerCentre.second - 2*radius; j < layerCentre.second + 2*radius; ++j){
                valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
                // check if the cell is in a sphere
                bool is_in_sphere = false;
                float x=0, y=0;
                convertFromGridIdxToCoordinates(x,y,i,j,wld);
                Eigen::Vector3f cellPos(x, y, wld->getAirLayer(k)->getHeight());
                for (unsigned int s_idx = 0; s_idx<close_spheres.size(); s_idx++)
                {
                    if ((cellPos - plume.free_spheres[close_spheres[s_idx]].center).norm() < plume.free_spheres[close_spheres[s_idx]].r)
                    {
                        is_in_sphere = true;
                        break;
                    }
                }
                if (dist <= 2*radius && is_in_sphere){
                    valType normDist = (valType) dist / ((valType) (2.*radius));
                    valType newTemp = PWM::Utils::interpolate(sLTemp, wld->getAirLayer(layerIndex)->getTemperature(i, j), normDist);
                    valType newVelX, newVelY;
                    if (airHeight < 8000){
                        newVelX = wld->getAirLayer(layerIndex)->getVelocityPhi(i, j);
                        newVelY = wld->getAirLayer(layerIndex)->getVelocityTheta(i, j);
                    }
                    else{
                        newVelX = wld->getAirLayer(layerIndex)->getVelocityPhi(i, j) * 0.01 * sLVertVel;
                        newVelY = wld->getAirLayer(layerIndex)->getVelocityTheta(i, j) * 0.01 * sLVertVel;
                    }
                    //valType newCW = PWM::Utils::interpolate<valType>(sLMoist, 0, normDist) + wld->getAirLayer(layerIndex)->getMoisture(i, j);
                    valType newCW = sLMoist + wld->getAirLayer(layerIndex)->getMoisture(i, j); // JG moisture
                    wld->getAirLayer(layerIndex)->setTemperature(i, j, newTemp);
                    wld->getAirLayer(layerIndex)->setVelocityPhi(i, j, newVelX);
                    wld->getAirLayer(layerIndex)->setVelocityTheta(i, j, newVelY);
                    wld->getAirLayer(layerIndex)->setMoisture(i, j, newCW);
                }
            }
        }
    }

    // old
    /*for (int k = 0; k < plume.smoke_layers.size(); ++k){
        int radius = std::round(plume.smoke_layers[k].r / wld->getCellSize());
        valType sLTemp = plume.smoke_layers[k].temperature;
        auto layerLoc = plume.smoke_layers[k].center;
        int layerIndex = findAirLayerIndex(wld->getAirLayers(), layerLoc.z());
        if (layerIndex == -1 || radius > (wld->getWorldSize() / 2))//-1 means that the smoke layer has passed above the highest atmosphere layer
            continue;
        // Easy enough to avoid by increasing the number of layers and atmotop, though the model behaviour becomes less accurate above 10 km
        auto layerCentre = convertFromCoordinatesToGridIdx((float) layerLoc.x(), (float) layerLoc.y(), wld);
        for (int i = layerCentre.first - radius; i < layerCentre.first + radius; ++i){
            for (int j = layerCentre.second - radius; j < layerCentre.first + radius; ++j){
                valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(layerCentre, std::make_pair(i, j));
                if (dist <= radius){
                    valType normDist = (valType) dist / (valType) radius;
                    valType newTemp = PWM::Utils::interpolate(sLTemp, wld->getAirLayer(layerIndex)->getTemperature(i, j), normDist);
                    valType newVelX = PWM::Utils::interpolate<valType>(0, wld->getAirLayer(layerIndex)->getVelocityPhi(i, j), normDist);
                    valType newVelY = PWM::Utils::interpolate<valType>(0, wld->getAirLayer(layerIndex)->getVelocityTheta(i, j), normDist);
                    valType newCW = PWM::Utils::interpolate<valType>(0, 1, normDist);
                    wld->getAirLayer(layerIndex)->setTemperature(i, j, newTemp);
                    wld->getAirLayer(layerIndex)->setVelocityPhi(i, j, newVelX);
                    wld->getAirLayer(layerIndex)->setVelocityTheta(i, j, newVelY);
                    wld->getAirLayer(layerIndex)->setMoisture(i, j, newCW);
                }
            }
        }
    }*/
}

//Testing the first validation case: donut plume.
int donutClouds(int steps, int width, int height, valType xSize, valType ySize, int layers, valType topOfAtmo){
    /*bool writeTestGrids = false;
    if (writeTestGrids){
        int width = 128;
        auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
        auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj", 400.f, 15000.f);
        y->maxLayers = 15;
        y->atmoTop = 18000.f;
        y->init("../resources/PlanetEarth.json", ter, width);
        for (int i = 0; i < y->getAirLayers().size(); ++i){
            y->getAirLayer(i)->getVelocityTheta() = 0;
            y->getAirLayer(i)->getVelocityPhi() = 0;
            y->getAirLayer(i)->getCondensedWater() = 0;
        }

        auto vdbExporter1 = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

        float weatherDT = 1.f;
        double plumeDT = 0.002;
        auto t = worldEngine(y, weatherDT, true);

        scene_model plume;
        plume.setup_data(ter, y->getPlanet(), t.getWorldModel()->getAirLayers().size(), plumeDT);
        // take altitude layers from atmospheric model
        for (unsigned int i = 0; i< t.getWorldModel()->getAirLayers().size(); i++)
        {
            plume.wind_altitudes[i] = t.getWorldModel()->getAirLayer(i)->getHeight();
        }
        plume.r_0 = 200.f;
        for (int i = 0; i < 20000; ++i){
            float d = plume.run_time;
            plume.frame_draw(plumeDT);
            transferDataFromPlumeToAtmAlt(plume, t.getWorldModel());
            if (i % ((int) (weatherDT / plumeDT)) == 0)
                t.step();
        }
        for (int i = 0; i < y->getAirLayers().size(); ++i)
            PWM::Utils::writeTerrElevImage("../output/PlumeLocatedDensity/Layer_" + std::to_string(i) + "_CloudDensity.ppm", y->getAirLayer(i)->getCondensedWater());
        vdbExporter1.write(t.getWorldModel(), "../output/PlumeLocatedDensity", 2000, t.getSimTimePassed());

        for (int i = 0; i < y->getAirLayers().size(); ++i){
            y->getAirLayer(i)->getVelocityPhi() = 1;//column - x-direction
            y->getAirLayer(i)->getVelocityTheta() = 0;//row - y-direction
            for (int j = 0; j < width; ++j){//row
                for (int k = 0; k < width; ++k){//column
                    valType val = (width - 1 - k) / (valType) width;
                    y->getAirLayer(i)->setCondensedWater(j, k, val);
                }
            }
        }
        for (int i = 0; i < y->getAirLayers().size(); ++i)
            PWM::Utils::writeTerrElevImage("../output/XDensityGradient/Layer_" + std::to_string(i) + "_CloudDensity.ppm", y->getAirLayer(i)->getCondensedWater());
        vdbExporter1.write(y, "../output/XDensityGradient", 0, 0);

        for (int i = 0; i < y->getAirLayers().size(); ++i){
            y->getAirLayer(i)->getVelocityPhi() = 0;//column - x-direction
            y->getAirLayer(i)->getVelocityTheta() = -1;//row - y-direction
            for (int j = 0; j < width; ++j){//row
                valType val = j / (valType) width;
                for (int k = 0; k < width; ++k){//column
                    y->getAirLayer(i)->setCondensedWater(j, k, val);
                }
            }
        }
        for (int i = 0; i < y->getAirLayers().size(); ++i)
            PWM::Utils::writeTerrElevImage("../output/YDensityGradient/Layer_" + std::to_string(i) + "_CloudDensity.ppm", y->getAirLayer(i)->getCondensedWater());
        vdbExporter1.write(y, "../output/YDensityGradient", 0, 0);
    }*/

    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/DonutCloud");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj", 400.f, xSize / 2);
//    auto ter = terrain_structure();
//    ter = ter.testTerrainStructure(128);

    valType cloudThreshold = 0.1;
    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;

    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, width, height, xSize, ySize, layerHeights);

    for (int i = 0; i < y->getAirLayers().size(); ++i){
        y->getAirLayer(i)->getVelocityTheta() = 0;
        y->getAirLayer(i)->getVelocityPhi() = 0;
        y->getAirLayer(i)->getCondensedWater() = 0;
        y->getAirLayer(i)->getMoisture() = 0.005;
//        y->getAirLayer(i)->getTemperature() = 300 - (100 * (10000 / y->getAirLayer(i)->getHeight()));
    }
    valType val = 0.105;

    int cloudIndex = y->getLayerAtHeight(3300);

    y->getAirLayer(cloudIndex)->getTemperature() = 245;
    y->getAirLayer(cloudIndex)->getCondensedWater() = val;

    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.ppm", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.ppm", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    float weatherDT = 2.f;
    float cveDT = 6.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(ter, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(0, 0);
    plume.setPlumeLoc(Eigen::Vector3f(0, 0, ventHeight));
//    float ventHeight = ter.sampleHeightAt(plumeLocation.first, plumeLocation.second);
//    plume.setPlumeLoc(Eigen::Vector3f(plumeLocation.first, plumeLocation.second, ventHeight));
//     take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(0, 0, y);
//    std::cout << "Plume z_0 is at " << ventHeight << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    //plume.r_0 = 200.f; <- NEVER DO THAT

    /* if (b.getWorldModel() == y)
        std::cout << "\033[1;32mUsing a real world map works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Using a real world map does not perform as expected!\033[0m" << std::endl;
     */
    double previousStep = 0;
    auto manCloudGen = [](const std::shared_ptr<PWM::Model::airLayer<dsType, valType>>& layer, const size_t gap, const valType val){
        #pragma omp parallel for
            for (size_t i = 0; i < layer->getObstacles().size(); ++i){
                if ((i / gap) % 5 == 0)
                    layer->setCondensedWater(i, val);
                else if ((i / gap) % 5 == 1)
                    layer->setCondensedWater(i, val / 100);
            }
        #pragma omp barrier
    };

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = 100;
    b.outputDirectory = "../ImageOutput/MergeTest/DonutCloud";
    for (int i = 0; i <= steps; ++i){
//		std::cerr << "Performing step " << i << " for the plume engine." << std::endl;
        float d = plume.run_time;
        plume.frame_draw(plumeDT);
//        transferDataFromPlumeToAtm(plume, b);
        d = plume.run_time - d;
        int frame_num = plume.frame_count / plume_export_frequency;
//		std::cout << "Step " << plume.frame_count << " of the plume engine took " << d << " seconds." << std::endl;
        if (i % plume_export_frequency == 0){
//            plume.export_spheres_binary(plume_export_frequency);
//            plume.export_spheres_combined(plume_export_frequency);
            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/DonutCloud", i, cloudThreshold);
//            vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold);
//            if (b.getStepCount() % 2  == 0){
                PWM::Utils::writeAirImages("../ImageOutput/MergeTest/DonutCloud", frame_num, b.getWorldModel());
//            }
        }
        if (i % ((int) (weatherDT / plumeDT)) == 0){
//            std::cerr << "Performing step " << b.getStepCount() << " for the weather engine." << std::endl;
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());
            previousStep = b.getRunTimePassed();
            b.step();
            transferDataFromAtmToPlume(plume, b);
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
//            if (i == 7800)
//                std::cerr << "debug point reached" << std::endl;
            if (i == steps / 8)
                std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
            else if (i == steps / 4)
                std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 8))
                std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
            else if (i == steps / 2)
                std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
            else if (i == (5 * steps / 8))
                std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 4))
                std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
            else if (i == (7 * steps / 8))
                std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
        }
    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cout << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cout << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

//Testing the second validation case: ash rain.
int ashRain(int steps, int width, int height, valType xSize, valType ySize, int layers, valType topOfAtmo, bool terse = false){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/AshRain");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj", 400.f, xSize / 2);
//    auto ter = terrain_structure();
//    ter = ter.testTerrainStructure(128);
    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    valType cloudThreshold = 0.1;
    y->init("../resources/PlanetEarth.json", ter, width, height, xSize, ySize, layerHeights);

    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed };

    const siv::PerlinNoise::seed_type moisSeed = 647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11883u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    int octaves = 8;

    for (int k = 0; k < y->getAirLayers().size(); ++k){
        if (y->getAirLayer(k)->getHeight() < 7000){
            y->getAirLayer(k)->getVelocityTheta() = 1. + (0.005 * y->getAirLayer(k)->getHeight());
            y->getAirLayer(k)->getVelocityPhi() = 0;
        }
        else{
            y->getAirLayer(k)->getVelocityTheta() = 0;
            y->getAirLayer(k)->getVelocityPhi() = 1. + (0.005 * y->getAirLayer(k)->getHeight());
        }
        y->getAirLayer(k)->getCondensedWater() = 0.05;
        y->getAirLayer(k)->getMoisture() = 0.01;
//        std::stringstream kTemp, kMois, kCloud;
//        kTemp << "Layer " << k << " temperature (Initial temp = " << y->getAirLayer(k)->getTemperature(0) << "):" << "\n";
//        kMois << "Layer " << k << " moisture:" << "\n";
//        kCloud << "Layer " << k << " clouds:" << "\n";
        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.05 : 0.12;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 4000) || (y->getAirLayer(k)->getHeight() > 10000) ? 0.02 : 0.08;
        for (int i = 0; i < height; ++i){
            for (int j = 0; j < width; ++j){
                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
//                kTemp << newTemp << "\t";
                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
//                kMois << newMois << "\t";
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
//                kCloud << newCloud << "\t";
            }
//            kTemp << "\n";
//            kMois << "\n";
//            kCloud << "\n";
        }
//        std::cout << kTemp.str() << "\n" << std::endl;
//        std::cout << kMois.str() << "\n" << std::endl;
//        std::cout << kCloud.str() << "\n" << std::endl;
    }

    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.ppm", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.ppm", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    float weatherDT = 2.f;
    float cveDT = 6.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(ter, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
//    float ventHeight = ter.sampleHeightAt(plumeLocation.first, plumeLocation.second);
//    plume.setPlumeLoc(Eigen::Vector3f(plumeLocation.first, plumeLocation.second, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    float ventHeight = ter.sampleHeightAt(0, 0);
    plume.setPlumeLoc(Eigen::Vector3f(0, 0, ventHeight));
//    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    //plume.r_0 = 200.f; <- NEVER DO THAT

    /* if (b.getWorldModel() == y)
        std::cout << "\033[1;32mUsing a real world map works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Using a real world map does not perform as expected!\033[0m" << std::endl;
     */
    double previousStep = 0;
    auto manCloudGen = [](const std::shared_ptr<PWM::Model::airLayer<dsType, valType>>& layer, const size_t gap, const valType val){
        #pragma omp parallel for
            for (size_t i = 0; i < layer->getObstacles().size(); ++i){
                if ((i / gap) % 5 == 0)
                    layer->setCondensedWater(i, val);
                else if ((i / gap) % 5 == 1)
                    layer->setCondensedWater(i, val / 100);
            }
        #pragma omp barrier
    };

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = 100;
    b.outputDirectory = "../ImageOutput/MergeTest/AshRain";
    for (int i = 0; i <= steps; ++i){
//		std::cerr << "Performing step " << i << " for the plume engine." << std::endl;
        float d = plume.run_time;
        plume.frame_draw(plumeDT);
//        transferDataFromPlumeToAtm(plume, b);
        d = plume.run_time - d;
        int frame_num = plume.frame_count / plume_export_frequency;
//		std::cout << "Step " << plume.frame_count << " of the plume engine took " << d << " seconds." << std::endl;
//        if (i % plume_export_frequency == 0){
//            plume.export_spheres_binary(plume_export_frequency, true);
            // plume.export_spheres_combined(plume_export_frequency);              
//        }
        if (i % ((int) (weatherDT / plumeDT)) == 0){
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());
            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/AshRain", i, cloudThreshold);
            if(!terse)
                vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold);
            if (b.getStepCount() % 2 == 0){
                if(terse)
                    PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/AshRain", i, b.getWorldModel());
                else
                    PWM::Utils::writeAirImages("../ImageOutput/MergeTest/AshRain", i, b.getWorldModel());
            }
//            std::cerr << "Performing step " << b.getStepCount() << " for the weather engine." << std::endl;
            previousStep = b.getRunTimePassed();
            b.step();
            transferDataFromAtmToPlume(plume, b);
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
            if (i == steps / 8)
                std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
            else if (i == steps / 4)
                std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 8))
                std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
            else if (i == steps / 2)
                std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
            else if (i == (5 * steps / 8))
                std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 4))
                std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
            else if (i == (7 * steps / 8))
                std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
        }
    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cerr << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cerr << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

//Testing the first validation case: donut plume, using the bigger scale landscape.
int donutCloudsBigScale(int steps, int layers, valType topOfAtmo){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/DonutCloudBigScale");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructBigTerrain("../resources/MtStHelens.elv");
    std::cout << "Terrain has " << ter.field_size << " cells a side, with a length of " << ter.max_xyz - ter.min_xyz << " metres on a side. Each cell is " << ter.cell_size << " metres on a side." << std::endl;

    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;

    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, layerHeights);

    for (int i = 0; i < y->getAirLayers().size(); ++i){
        y->getAirLayer(i)->getVelocityTheta() = 0;
        y->getAirLayer(i)->getVelocityPhi() = 0;
        y->getAirLayer(i)->getCondensedWater() = 0;
        y->getAirLayer(i)->getMoisture() = 0.005;
        //        y->getAirLayer(i)->getTemperature() = 300 - (100 * (10000 / y->getAirLayer(i)->getHeight()));
    }
    valType val = 0.105;

    valType cloudThreshold = 0.1;
    int cloudIndex = y->getLayerAtHeight(3300);

    y->getAirLayer(cloudIndex)->getTemperature() = 245;
    y->getAirLayer(cloudIndex)->getCondensedWater() = val;

    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.ppm", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.ppm", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    float weatherDT = 2.f;
    float cveDT = 6.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(ter, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);
    plume.setPlumeLoc(Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    //plume.r_0 = 200.f; <- NEVER DO THAT

    /* if (b.getWorldModel() == y)
        std::cout << "\033[1;32mUsing a real world map works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Using a real world map does not perform as expected!\033[0m" << std::endl;
     */
    double previousStep = 0;
    auto manCloudGen = [](const std::shared_ptr<PWM::Model::airLayer<dsType, valType>>& layer, const size_t gap, const valType val){
        #pragma omp parallel for
        for (size_t i = 0; i < layer->getObstacles().size(); ++i){
            if ((i / gap) % 5 == 0)
                layer->setCondensedWater(i, val);
            else if ((i / gap) % 5 == 1)
                layer->setCondensedWater(i, val / 100);
        }
        #pragma omp barrier
    };

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = 100;
    b.outputDirectory = "../ImageOutput/MergeTest/DonutCloudBigScale";
    for (int i = 0; i <= steps; ++i){
        //		std::cerr << "Performing step " << i << " for the plume engine." << std::endl;
        float d = plume.run_time;
        plume.frame_draw(plumeDT);
        //        transferDataFromPlumeToAtm(plume, b);
        transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
        transferAshFromPlumeToAtm(plume, b.getWorldModel());
        d = plume.run_time - d;
        int frame_num = plume.frame_count / plume_export_frequency;
        //		std::cout << "Step " << plume.frame_count << " of the plume engine took " << d << " seconds." << std::endl;
        if (i % plume_export_frequency == 0){
            // plume.export_spheres(plume_export_frequency);
            plume.export_spheres_binary(plume_export_frequency);
            // plume.export_spheres_combined(plume_export_frequency);
            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/DonutCloud", i, cloudThreshold);
            vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold);

            if (b.getStepCount() % 2  == 0){
                PWM::Utils::writeAirImages("../ImageOutput/MergeTest/DonutCloudBigScale", b.getStepCount(), b.getWorldModel());
            }
        }
        if (i % ((int) (weatherDT / plumeDT)) == 0){
            //            std::cerr << "Performing step " << b.getStepCount() << " for the weather engine." << std::endl;
            previousStep = b.getRunTimePassed();
            b.step();
            transferDataFromAtmToPlume(plume, b);
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
            if (i == 7800)
                std::cerr << "debug point reached" << std::endl;
            if (i == steps / 8)
                std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
            else if (i == steps / 4)
                std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 8))
                std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
            else if (i == steps / 2)
                std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
            else if (i == (5 * steps / 8))
                std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 4))
                std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
            else if (i == (7 * steps / 8))
                std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
        }
    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cout << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cout << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

//Testing the second validation case: ash rain, with the bigger landscape.
int ashRainBigScale(int steps, int layers, valType topOfAtmo, bool terse = false, bool analysis = false, bool atmosim = true){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/AshRainBigScale");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructBigTerrain("../resources/MtStHelens.elv");
    std::cout << "Terrain has " << ter.field_size << " cells a side, with a length of " << ter.max_xyz - ter.min_xyz << " metres on a side. Each cell is " << ter.cell_size << " metres on a side." << std::endl;

    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, layerHeights);
    
    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed };

    const siv::PerlinNoise::seed_type moisSeed = 647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11883u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    int octaves = 8;

    /*
    std::cerr << "TEMPERATURES" << std::endl;
    for (int k = 0; k < y->getAirLayers().size(); ++k)
        std::cerr << "layer " << k << " temp = " << y->getAirLayer(k)->getTemperature(512, 512) << std::endl;
     */

    
    std::cout << "Adding noise to initial values." << std::endl;
    for (int k = 0; k < y->getAirLayers().size(); ++k){
        if (y->getAirLayer(k)->getHeight() < 7000){
            y->getAirLayer(k)->getVelocityTheta() = 1. + (0.005 * y->getAirLayer(k)->getHeight());
            y->getAirLayer(k)->getVelocityPhi() = 0;
        }
        else{
            y->getAirLayer(k)->getVelocityTheta() = 0;
            y->getAirLayer(k)->getVelocityPhi() = 1. + (0.005 * y->getAirLayer(k)->getHeight());
        }
        //        std::stringstream kTemp, kMois, kCloud;
        //        kTemp << "Layer " << k << " temperature (Initial temp = " << y->getAirLayer(k)->getTemperature(0) << "):" << "\n";
        //        kMois << "Layer " << k << " moisture:" << "\n";
        //        kCloud << "Layer " << k << " clouds:" << "\n";
        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.05 : 0.12;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 4000) || (y->getAirLayer(k)->getHeight() > 10000) ? 0.02 : 0.08;
        for (int i = 0; i < ter.field_size; ++i){
            for (int j = 0; j < ter.field_size; ++j){
                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
                //                kTemp << newTemp << "\t";
                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
                //                kMois << newMois << "\t";
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
                //                kCloud << newCloud << "\t";
            }
            //            kTemp << "\n";
            //            kMois << "\n";
            //            kCloud << "\n";
        }
//                std::cout << kTemp.str() << std::endl;
        //        std::cout << kMois.str() << "\n" << std::endl;
        //        std::cout << kCloud.str() << "\n" << std::endl;
    }
    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.ppm", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.ppm", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    valType cloudThreshold = 0.1;
    
    float weatherDT = 4.f;
    float cveDT = 12.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(ter, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);
    
    plume.setPlumeLoc(Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    std::cout << "Plume r_0 is " << plume.r_0 << " metre." << std::endl;
    //plume.r_0 = 200.f; <- NEVER DO THAT

    /* if (b.getWorldModel() == y)
        std::cout << "\033[1;32mUsing a real world map works as expected.\033[0m" << std::endl;
    else
        std::cout << "\033[1;31mError! Using a real world map does not perform as expected!\033[0m" << std::endl;
     */
    double previousStep = 0, d = 0;
    auto manCloudGen = [](const std::shared_ptr<PWM::Model::airLayer<dsType, valType>>& layer, const size_t gap, const valType val){
        #pragma omp parallel for
        for (size_t i = 0; i < layer->getObstacles().size(); ++i){
            if ((i / gap) % 5 == 0)
                layer->setCondensedWater(i, val);
            else if ((i / gap) % 5 == 1)
                layer->setCondensedWater(i, val / 100);
        }
        #pragma omp barrier
    };

     // Fake hotspot
    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
    int plumeWidth =  2 * std::ceil(plume.r_0 / ter.cell_size);
    auto pl = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz - plumeWidth, BigPlumeLocation.second + ter.min_xyz - plumeWidth, y);
    int startj = pl.first, startk = pl.second;
    double highT = plume.T_0;
    for (int j = -plumeWidth; j < plumeWidth; ++j){
        for (int k = -plumeWidth; k < plumeWidth; ++k){
            auto cellLoc = std::make_pair(j + startj, k + startk);
            valType dist = PWM::Utils::calcCartesianDistance<valType>(plumeLoc, cellLoc);
            if (dist <= plumeWidth){
                y->getTerrain()->setTemperature(ijToI(ter.field_size, j + startj, k + startk), highT);
            }
        }
    }

    PWM::Utils::writeTempImage("../ImageOutput/MergeTest/terrainModelTemp.ppm", y->getTerrain()->getTemperature());

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = 100;
    b.outputDirectory = "../ImageOutput/MergeTest/AshRainBigScale";
    for (int i = 0; i <= steps; ++i){
        //		std::cerr << "Performing step " << i << " for the plume engine." << std::endl;

        float d = plume.run_time;
        Timer t;
        
        // t.start();

        float currD = plume.run_time;

        plume.frame_draw(plumeDT);
        // t.stop();
        // std::cout << "Plume " << i << " took " << t.peek() << " s" << std::endl;
      
      
        //        transferDataFromPlumeToAtm(plume, b);
        int frame_num = plume.frame_count / plume_export_frequency;
        //		std::cout << "Step " << plume.frame_count << " of the plume engine took " << d << " seconds." << std::endl;
        if (i % plume_export_frequency == 0){
            if(!analysis)
            {
                // t.start();
                plume.export_spheres_binary(plume_export_frequency, true);
             
                // t.stop();
                // std::cout << "Skirt & plume output took " << t.peek() << "s" << std::endl;
                // plume.export_spheres_combined(plume_export_frequency);
            }
        }
        if (i % ((int) (weatherDT / plumeDT)) == 0){
            t.start();
            std::cout << "WEATHER PHASE" << std::endl;
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());
            std::cout << "DATA TRANSFER" << std::endl;
            if(atmosim)
            {
                PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/AshRainBigScale", i, cloudThreshold);
                if(!terse && !analysis)
                    vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold);
                if (b.getStepCount() % 2 == 0){
                    if(!analysis)
                    {
                        if(terse)
                            PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/AshRainBigScale", i, b.getWorldModel());
                        else
                            PWM::Utils::writeAirImages("../ImageOutput/MergeTest/AshRainBigScale", i, b.getWorldModel());
                       
                    }
                }
            }
           
            // std::cout << "SKIRT OUTPUT" << std::endl;
            //            std::cerr << "Performing step " << b.getStepCount() << " for the weather engine." << std::endl;
            // previousStep = b.getRunTimePassed();
            if(atmosim)
                b.step();
            std::cout << "ATMOSIM" << std::endl;
            
            t.stop();
            std::cout << "Weather phase took " << t.peek() << "s" << std::endl;
            
            // transferDataFromAtmToPlume(plume, b);
            // double currentStep = b.getRunTimePassed();
            // std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;

            if (i == steps / 8)
                std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
            else if (i == steps / 4)
                std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 8))
                std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
            else if (i == steps / 2)
                std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
            else if (i == (5 * steps / 8))
                std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 4))
                std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
            else if (i == (7 * steps / 8))
                std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
        }
    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cout << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cout << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

//Testing the second validation case: ash rain, with the bigger landscape.
int ashRainBigScaleFinal(int steps, int layers, valType topOfAtmo, int startSave = 0){
    Timer t;
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/AshRainBigScale");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".txt" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructHalfTerrain("../resources/MtStHelens.elv");
    auto pTer = constructHalfTerrain("../resources/MtStHelens.elv");
    std::cout << "Terrain has " << ter.field_size << " cells a side, with a length of " << ter.max_xyz - ter.min_xyz << " metres on a side. Each cell is " << ter.cell_size << " metres on a side." << std::endl;

    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    valType cloudThreshold = 0.1;
    y->init("../resources/PlanetEarth.json", ter, layerHeights);
    
    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed};

    const siv::PerlinNoise::seed_type moisSeed = 647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11883u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    const siv::PerlinNoise::seed_type velXSeed = 362416468u;
    const siv::PerlinNoise velXPerNoiGen{ velXSeed };

    const siv::PerlinNoise::seed_type velYSeed = 1231578445u;
    const siv::PerlinNoise velYPerNoiGen{ velYSeed };

    int octaves = 8;

    std::cout << "Adding noise to initial values." << std::endl;
    for (int k = 0; k < y->getAirLayers().size(); ++k){
        valType velXBase, velYBase;
        // if (y->getAirLayer(k)->getHeight() < 7000){
            // velXBase = 10;
            // velYBase = 0;
        // }
        // else{
        // y->getAirLayer(k)->getVelocityPhi() = velXBase = 10;//. + (0.001 * y->getAirLayer(k)->getHeight());
        // y->getAirLayer(k)->getVelocityTheta() = velYBase = 0;
        // }

        valType R = 6371000, T = 86400 * 12, t = 1;
        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.08 : 0.06;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 4000) || (y->getAirLayer(k)->getHeight() > 10000) ? 0.06 : 0.08;
        auto interval = std::numbers::pi / ter.field_size;
        for (int i = 0; i < ter.field_size; ++i){
            for (int j = 0; j < ter.field_size; ++j){
                valType thetaCoor, lambdaPrime;
                thetaCoor = i * interval;
                lambdaPrime = j * interval - (2 * std::numbers::pi * t) / T;
                valType velPhi = -(10 * R / T) * std::sin(2 * lambdaPrime) * std::cos(thetaCoor) * std::cos((std::numbers::pi * t) / T);
                y->getAirLayer(k)->setVelocityPhi(i, j, velPhi * 10);
                valType velTheta = -(10 * R / T) * std::pow(std::sin(lambdaPrime), 2) * std::sin(2 * thetaCoor) * std::cos((std::numbers::pi * t) / T) + (2 * std::numbers::pi * R / T) * std::cos(thetaCoor);
                y->getAirLayer(k)->setVelocityTheta(i, j, velTheta * 10);

                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
                //                kTemp << newTemp << "\t";

                // valType newVelX = velXBase;// + (velXPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 10);
                // y->getAirLayer(k)->setVelocityPhi(i, j, newVelX);

                // valType newVelY = velYBase;// + (velYPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 10);
                // y->getAirLayer(k)->setVelocityTheta(i, j, newVelY);

                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
                //                kMois << newMois << "\t";
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
                //                kCloud << newCloud << "\t";
            }
        }
    }
    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.png", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.png", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();
    writeLayerStruct("../ImageOutput/MergeTest/AshRainBigScale", y);

    float weatherDT = 3.f;
    float cveDT = 9.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(pTer, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);
    writeLayerStruct("../ImageOutput/MergeTest/AshRainBigScale", y);

    plume.setPlumeLoc(Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    std::cout << "Plume r_0 is " << plume.r_0 << " metre." << std::endl;
   
    double previousStep = 0;

//    // Fake hotspot
//    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
//    int plumeWidth =  2 * std::ceil(plume.r_0 / ter.cell_size);
//    auto pl = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz - plumeWidth, BigPlumeLocation.second + ter.min_xyz - plumeWidth, y);
//    int startj = pl.first, startk = pl.second;
//    double highT = plume.T_0 / 2;
//    for (int j = -plumeWidth; j < plumeWidth; ++j){
//        for (int k = -plumeWidth; k < plumeWidth; ++k){
//            auto cellLoc = std::make_pair(j + startj, k + startk);
//            valType dist = PWM::Utils::calcCartesianDistance<valType>(plumeLoc, cellLoc);
//            if (dist <= plumeWidth){
//                y->getTerrain()->setTemperature(ijToI(ter.field_size, j + startj, k + startk), highT);
//            }
//        }
//    }

    PWM::Utils::writeTempImage("../ImageOutput/MergeTest/terrainModelTemp.ppm", y->getTerrain()->getTemperature());

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = weatherDT / plumeDT;
    b.outputDirectory = "../ImageOutput/MergeTest/AshRainBigScale";
    for (int i = 0; i <= steps; ++i){
        t.start();
        float d = plume.run_time;
        plume.frame_draw(plumeDT);
        //        transferDataFromPlumeToAtm(plume, b);
        d = plume.run_time - d;
        int frame_num = plume.frame_count / plume_export_frequency;
        
        if (i % plume_export_frequency == 0 && i >= startSave * plume_export_frequency){
            plume.export_spheres_binary(plume_export_frequency, true);
        }
        if (i % ((int) (weatherDT / plumeDT)) == 0){
            PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/AshRainBigScale", i, b.getWorldModel());
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());
            if(i >= startSave * plume_export_frequency){
               PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/AshRainBigScale", i, cloudThreshold);
               // vdbExporter.write(b.getWorldModel(), "../output", b.getStepCount(), b.getSimTimePassed(), cloudThreshold);
               writeTerrainAshVals("../ImageOutput/MergeTest/AshRainBigScale", b.getWorldModel()->getTerrain(), i, plumeLoc);
            }
            previousStep = b.getRunTimePassed();
            b.step();
            transferDataFromAtmToPlume(plume, b);
            stripEdgeAsh(b.getWorldModel()->getAirLayers(), b.getWorldModel()->getConvectionLayers());
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
        }
        t.stop();
        std::cerr << "Step " << i << " took " << t.peek() << "s" << std::endl;
        if (i == steps / 8)
            std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
        else if (i == steps / 4)
            std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 8))
            std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
        else if (i == steps / 2)
            std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
        else if (i == (5 * steps / 8))
            std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 4))
            std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
        else if (i == (7 * steps / 8))
            std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;

    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cout << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cout << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

//Testing the ash validation case: ash rain on terrain, with much bigger landscape.
int ashRainValidation(int steps, int layers, valType topOfAtmo, int cellCount, float landScapeMax, int startSave = 0){
    Timer t;
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/AshRainValidation");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".txt" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructFlatBigTerrain(landScapeMax, cellCount);
    auto pTer = constructFlatBigTerrain(landScapeMax, cellCount);
    std::cout << "Terrain has " << ter.field_size << " cells a side, with a length of " << ter.max_xyz - ter.min_xyz << " metres on a side. Each cell is " << ter.cell_size << " metres on a side." << std::endl;

    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    valType cloudThreshold = 0.1;
    y->init("../resources/PlanetEarth.json", ter, layerHeights);

    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed};

    const siv::PerlinNoise::seed_type moisSeed = 647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11883u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    const siv::PerlinNoise::seed_type velXSeed = 362416468u;
    const siv::PerlinNoise velXPerNoiGen{ velXSeed };

    const siv::PerlinNoise::seed_type velYSeed = 1231578445u;
    const siv::PerlinNoise velYPerNoiGen{ velYSeed };

    int octaves = 8;

    std::cout << "Adding noise to initial values." << std::endl;
    for (int k = 0; k < y->getAirLayers().size(); ++k){
        valType velXBase, velYBase;
        if (y->getAirLayer(k)->getHeight() < 7000){
            velXBase = 0;
            velYBase = 1. + (0.001 * y->getAirLayer(k)->getHeight());
        }
        else{
            velXBase = 1. + (0.001 * y->getAirLayer(k)->getHeight());
            velYBase = 0;
        }

        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.08 : 0.06;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 4000) || (y->getAirLayer(k)->getHeight() > 10000) ? 0.06 : 0.08;
        for (int i = 0; i < ter.field_size; ++i){
            for (int j = 0; j < ter.field_size; ++j){
                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
                //                kTemp << newTemp << "\t";

                valType newVelX = velXBase;// + (velXPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 10);
                y->getAirLayer(k)->setVelocityPhi(i, j, newVelX);

                valType newVelY = velYBase;// + (velYPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 10);
                y->getAirLayer(k)->setVelocityTheta(i, j, newVelY);

                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
                //                kMois << newMois << "\t";
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
                //                kCloud << newCloud << "\t";
            }
        }
    }
    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.png", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.png", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();
    writeLayerStruct("../ImageOutput/MergeTest/AshRainValidation", y);

    float weatherDT = 3.f;
    float cveDT = 9.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(pTer, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);
    writeLayerStruct("../ImageOutput/MergeTest/AshRainValidation", y);

    plume.setPlumeLoc(Eigen::Vector3f(0, 0, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(0, 0, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    std::cout << "Plume r_0 is " << plume.r_0 << " metre." << std::endl;

    double previousStep = 0;

    //    // Fake hotspot
    //    auto ijToI = [](int width, int i, int j)->int{return i * width + j;};
    //    int plumeWidth =  2 * std::ceil(plume.r_0 / ter.cell_size);
    //    auto pl = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz - plumeWidth, BigPlumeLocation.second + ter.min_xyz - plumeWidth, y);
    //    int startj = pl.first, startk = pl.second;
    //    double highT = plume.T_0 / 2;
    //    for (int j = -plumeWidth; j < plumeWidth; ++j){
    //        for (int k = -plumeWidth; k < plumeWidth; ++k){
    //            auto cellLoc = std::make_pair(j + startj, k + startk);
    //            valType dist = PWM::Utils::calcCartesianDistance<valType>(plumeLoc, cellLoc);
    //            if (dist <= plumeWidth){
    //                y->getTerrain()->setTemperature(ijToI(ter.field_size, j + startj, k + startk), highT);
    //            }
    //        }
    //    }

    PWM::Utils::writeTempImage("../ImageOutput/MergeTest/terrainModelTemp.png", y->getTerrain()->getTemperature());

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = weatherDT / plumeDT;
    b.outputDirectory = "../ImageOutput/MergeTest/AshRainValidation";
    for (int i = 0; i <= steps; ++i){
        t.start();
        float d = plume.run_time;
        plume.frame_draw(plumeDT);
        //        transferDataFromPlumeToAtm(plume, b);
        d = plume.run_time - d;
        int frame_num = plume.frame_count / plume_export_frequency;

        if (i % plume_export_frequency == 0 && i >= startSave * plume_export_frequency){
            // plume.export_spheres_binary(plume_export_frequency, true);
        }
        if (i % ((int) (weatherDT / plumeDT)) == 0){
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());
            if(i >= startSave * plume_export_frequency){
                //                PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/AshRainValidation", i, cloudThreshold);
                //                vdbExporter.write(b.getWorldModel(), "../output", b.getStepCount(), b.getSimTimePassed(), cloudThreshold);
                // PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/AshRainValidation", i, b.getWorldModel());
                writeTerrainAshVals("../ImageOutput/MergeTest/AshRainValidation", b.getWorldModel()->getTerrain(), i, plumeLoc);
            }
            previousStep = b.getRunTimePassed();
            b.step();
            transferDataFromAtmToPlume(plume, b);
            stripEdgeAsh(b.getWorldModel()->getAirLayers(), b.getWorldModel()->getConvectionLayers());
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
        }
        t.stop();
        std::cerr << "Step " << i << " took " << t.peek() << "s" << std::endl;
        if (i == steps / 8)
            std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
        else if (i == steps / 4)
            std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 8))
            std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
        else if (i == steps / 2)
            std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
        else if (i == (5 * steps / 8))
            std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 4))
            std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
        else if (i == (7 * steps / 8))
            std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;

    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cout << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cout << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

// Shorter term simulation that includes three layers of skirts
int skirts3(int steps, int layers, valType topOfAtmo){
    Timer t;
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/Skirt");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructHalfTerrain("../resources/MtStHelens.elv");
   
    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;
    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, layerHeights);
    
    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed };

    const siv::PerlinNoise::seed_type moisSeed = 647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11883u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    int octaves = 8;

    std::cout << "Adding noise to initial values." << std::endl;
    for (int k = 0; k < y->getAirLayers().size(); ++k){
        if (y->getAirLayer(k)->getHeight() < 7000){
            y->getAirLayer(k)->getVelocityTheta() = 1. + (0.002 * y->getAirLayer(k)->getHeight());
            y->getAirLayer(k)->getVelocityPhi() = 0;
        }
        else{
            y->getAirLayer(k)->getVelocityTheta() = 0;
            y->getAirLayer(k)->getVelocityPhi() = 1. + (0.002 * y->getAirLayer(k)->getHeight());
        }
        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.075 : 0.1;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 6000) ? 0.16 : 0.08;
        for (int i = 0; i < ter.field_size; ++i){
            for (int j = 0; j < ter.field_size; ++j){
                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
            }
        }
    }
    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.ppm", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.ppm", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();

    valType cloudThreshold = 0.1;
    float weatherDT = 2.f;
    float cveDT = 6.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(ter, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);
    
    plume.setPlumeLoc(Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    std::cout << "Plume r_0 is " << plume.r_0 << " metre." << std::endl;
    
    // create skirt structures
    std::vector<skirt *> skirts;
    int numskirts = 3;
    Eigen::Vector3f center = Eigen::Vector3f(0.0f, -12000.0f, 5000.0f);
    
    for(int s = 0; s < numskirts; s++)
    {
        skirt * sk = new skirt();
        
        sk->init(center, 6500.0f, 150.0f, 11000.0f+(750.0f*(float)s), 8000.0f+(500.0f*(float)s), 100.0f, y);
        // JG PARAM radius = 6500, step between samples = 150, target height = 11000+layer diff, condensation height = 8000+layer diff, thickness = 100
        skirts.push_back(sk);
        center.z() += 1000.0f;
    }

    PWM::Utils::writeTempImage("../ImageOutput/MergeTest/terrainModelTemp.png", y->getTerrain()->getTemperature());

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = weatherDT / plumeDT;
    b.outputDirectory = "../ImageOutput/MergeTest/Skirt";
    for (int i = 0; i <= steps; ++i){
        t.start();
        
        plume.frame_draw(plumeDT);
      
        for(int s = 0; s < numskirts; s++)
        {
            skirts[s]->uplift(plume, plumeDT, 4000.0f, 2000.0f); // JG PARAM - horizontal influence zone = 4000, vertical influence zone = 2000
            skirts[s]->phaseTransition();
        }

        int frame_num = plume.frame_count / plume_export_frequency;

        if (i % plume_export_frequency == 0)
        {
            plume.export_spheres_binary(plume_export_frequency, true);
            for(int s = 0; s < numskirts; s++)
                skirts[s]->exportSkirt(frame_num, s);
        }

        if (i % ((int) (weatherDT / plumeDT)) == 0){
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());

            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/Skirt", i, cloudThreshold);
            vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold); // uncomment for final run
            PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/Skirt", i, b.getWorldModel()); // images of atmo layers

            b.step();
            transferDataFromAtmToPlume(plume, b);
        }
        t.stop();
        std::cerr << "Step " << i << " took " << t.peek() << "s" << std::endl;
        if (i == steps / 8)
            std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
        else if (i == steps / 4)
            std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 8))
            std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
        else if (i == steps / 2)
            std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
        else if (i == (5 * steps / 8))
            std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 4))
            std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
        else if (i == (7 * steps / 8))
            std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
    }

    b.printTiming(true);
    for(int s = 0; s < numskirts; s++)
        delete skirts[s];
    return 0;
}

// Shorter term simulation that includes a cap cloud and a Wilson cloud
int capWilson(int steps, int layers, valType topOfAtmo){
    Timer t;
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/CapWilson");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".txt" || dir_Entry.path().extension().string() == ".vdb"))
            std::filesystem::remove(dir_Entry.path());
    }

    currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("output/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
            std::filesystem::remove(dir_Entry.path());
    }
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructPartTerrain("../resources/MtStHelens.elv", 2);
    auto pTer = constructBigTerrain("../resources/MtStHelens.elv");

    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;
    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, layerHeights);

    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed };

    const siv::PerlinNoise::seed_type moisSeed = 1647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11898743u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    int octaves = 8;

    std::cout << "Adding noise to initial values." << std::endl;
    for (int k = 0; k < y->getAirLayers().size(); ++k){
        if (y->getAirLayer(k)->getHeight() < 7000){
            y->getAirLayer(k)->getVelocityTheta() = 1. + (0.002 * y->getAirLayer(k)->getHeight());
            y->getAirLayer(k)->getVelocityPhi() = 0;
        }
        else{
            y->getAirLayer(k)->getVelocityTheta() = 0;
            y->getAirLayer(k)->getVelocityPhi() = 1. + (0.002 * y->getAirLayer(k)->getHeight());
        }
        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.06 : 0.06;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 6000) ? 0.09 : 0.08;
        for (int i = 0; i < ter.field_size; ++i){
            for (int j = 0; j < ter.field_size; ++j){
                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
            }
        }
    }
    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.png", ter);
    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.png", y->getTerrain()->getElevation());
    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();
    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);

    // Create explosion struct
    auto explosion = std::make_shared<PWM::Model::pressureWave<valType>>(PWM::Model::pressureWave<valType>(y->getPlanet(), 1, 100000));
    explosion->setMaxLife(5);
    explosion->setCentre(BigPlumeLocation);
    explosion->setHeight(ventHeight + 500);

    valType cloudThreshold = 0.1;
    float weatherDT = 0.5f;
    float cveDT = 0.5f;
    float quickStepDT = 0.1f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);

    scene_model plume;
    plume.setup_data(pTer, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);

    plume.setPlumeLoc(Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
    std::cout << "Plume r_0 is " << plume.r_0 << " metre." << std::endl;

    // create skirt structures
    std::vector<skirt *> skirts;
    int numskirts = 3;
    Eigen::Vector3f center = Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, 9000.0f);

    for(int s = 0; s < numskirts; s++)
    {
        skirt * sk = new skirt();

        sk->init(center, 5000.0f, 200.0f, 15000.0f, 12000.0f, 100.0f, y);
        // JG PARAM radius = 8000, step between samples = 200, target height = 15000, condensation height = 12000, thickness = 100
        skirts.push_back(sk);
        center.z() += 200.0f; // hardly any seperation
    }
    center = Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight + 800);
    wilsonCloud* wc = new wilsonCloud(center, 2500.0, 50.0, 1, 5, y);
    wc->addExplosion(explosion);
    skirts.push_back((skirt*) wc);

    double previousStep = 0, d = 0;
    PWM::Utils::writeTempImage("../ImageOutput/MergeTest/terrainModelTemp.png", y->getTerrain()->getTemperature());

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = weatherDT / plumeDT;
    b.outputDirectory = "../ImageOutput/MergeTest/CapWilson";
    for (int i = 0; i <= steps; ++i){
        t.start();

        plume.frame_draw(plumeDT);

        for(int s = 0; s < skirts.size(); s++)
        {
            skirts[s]->uplift(plume, plumeDT, 8000.0f, 2000.0f); // JG PARAM - very broad influence zone horizontally, no tail off
            skirts[s]->phaseTransition();
        }
        wc->updateExp(plumeDT);

        int frame_num = plume.frame_count / plume_export_frequency;

        if (i % ((int) (weatherDT / plumeDT)) == 0){
            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/CapWilson", i, cloudThreshold);
            // PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/CapWilson", i, b.getWorldModel()); // images of atmo layers
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());

            previousStep = b.getRunTimePassed();
            b.step();
            stripEdgeAsh(b.getWorldModel()->getAirLayers(), b.getWorldModel()->getConvectionLayers());
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
            transferDataFromAtmToPlume(plume, b);
            for (auto s : skirts)
                s->updateTemp();
        }
        if (i % plume_export_frequency == 0)
        {
            plume.export_spheres_binary(plume_export_frequency, true);
            for(int s = 0; s < skirts.size(); s++)
                skirts[s]->exportSkirt(frame_num, s);

            vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold); // uncomment for final run
        }

        t.stop();
        std::cerr << "Step " << i << " took " << t.peek() << "s" << std::endl;
        if (i == steps / 8)
            std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
        else if (i == steps / 4)
            std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 8))
            std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
        else if (i == steps / 2)
            std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
        else if (i == (5 * steps / 8))
            std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
        else if (i == (3 * steps / 4))
            std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
        else if (i == (7 * steps / 8))
            std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
    }

    b.printTiming(true);
    for(int s = 0; s < numskirts; s++)
        delete skirts[s];
    return 0;
}

/*//auto currDir = std::filesystem::current_path();
//currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/CapWilson");
//for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
//    if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".txt" || dir_Entry.path().extension().string() == ".vdb"))
//        std::filesystem::remove(dir_Entry.path());
//}

//// Shorter term simulation that includes a cap cloud and Wilson clouds
//int capWilson(int steps, int layers, valType topOfAtmo){
//    Timer t;
//    auto currDir = std::filesystem::current_path();
//    currDir.remove_filename().remove_filename().concat("ImageOutput/").concat("MergeTest/CapWilson");
//    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
//        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".ppm" || dir_Entry.path().extension().string() == ".png" || dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".txt" || dir_Entry.path().extension().string() == ".vdb"))
//            std::filesystem::remove(dir_Entry.path());
//    }

//    currDir = std::filesystem::current_path();
//    currDir.remove_filename().remove_filename().concat("output/");
//    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
//        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && (dir_Entry.path().extension().string() == ".dat" || dir_Entry.path().extension().string() == ".vdb" || dir_Entry.path().extension().string() == ".txt"))
//            std::filesystem::remove(dir_Entry.path());
//    }
//    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
//    auto ter = constructPartTerrain("../resources/MtStHelens.elv", 2);
//    auto pTer = constructBigTerrain("../resources/MtStHelens.elv");
   
//    y->maxLayers = layers;
//    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
//    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

//    valType cloudThreshold = 0.1;

//    y->init("../resources/PlanetEarth.json", ter, layerHeights);
    
//    const siv::PerlinNoise::seed_type tempSeed = 100u;
//    const siv::PerlinNoise temPerNoiGen{ tempSeed };

//    const siv::PerlinNoise::seed_type moisSeed = 645527u;
//    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

//    const siv::PerlinNoise::seed_type cloudSeed = 186543u;
//    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

//    int octaves = 8;

//    std::cout << "Adding noise to initial values." << std::endl;
//    for (int k = 0; k < y->getAirLayers().size(); ++k){
//        y->getAirLayer(k)->getVelocityTheta() = 0.001 * y->getAirLayer(k)->getHeight();
//        y->getAirLayer(k)->getVelocityPhi() = 0;

//        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.075 : 0.1;
//        auto cloudVal = (y->getAirLayer(k)->getHeight() < 6000) ? 0.16 : 0.08;
//        for (int i = 0; i < ter.field_size; ++i){
//            for (int j = 0; j < ter.field_size; ++j){
//                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
//                y->getAirLayer(k)->setTemperature(i, j, newTemp);
//                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
//                y->getAirLayer(k)->setMoisture(i, j, newMois);
//                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
//                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
//            }
//        }
//    }
//    PWM::Utils::writeTerrElevImage<double>("../ImageOutput/MergeTest/terrainStructElev.ppm", ter);
//    PWM::Utils::writeTerrElevImage("../ImageOutput/MergeTest/terrainModelElev.ppm", y->getTerrain()->getElevation());
//    auto vdbExporter = PWM::Utils::vdbExporter<dsType, dsSType, valType, valType2>();
//    writeLayerStruct("../ImageOutput/MergeTest/CapWilson", y);

//    float ventHeight = ter.sampleHeightAt(BigPlumeLocation.first, BigPlumeLocation.second);
//    // Create explosion struct
//    auto explosion = std::make_shared<PWM::Model::pressureWave<valType>>(PWM::Model::pressureWave<valType>(y->getPlanet(), 1, 100000));
//    explosion->setMaxLife(15);
//    explosion->setCentre(BigPlumeLocation);
//    explosion->setHeight(ventHeight);
//    y->setExplosion(explosion);

//    float weatherDT = 2.0f;
//    float cveDT = 6.0f;
//    float quickStepDT = 0.1f;
//    double plumeDT = 0.02;
//    auto b = worldEngine(y, weatherDT, cveDT, true);

//    scene_model plume;
//    plume.setup_data(pTer, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
//    plume.addAtmoModel(y);
    
//    plume.setPlumeLoc(Eigen::Vector3f(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, ventHeight));
//    // take altitude layers from atmospheric model
//    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
//    {
//        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
//    }
//    auto plumeLoc = convertFromCoordinatesToGridIdx(BigPlumeLocation.first + ter.min_xyz, BigPlumeLocation.second + ter.min_xyz, y);
//    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;
//    std::cout << "Plume r_0 is " << plume.r_0 << " metre." << std::endl;
    
//    // create skirt structures for the cap cloud
//    std::vector<skirt *> skirts;
//    int numskirts = 3;
//    Eigen::Vector3f center = Eigen::Vector3f(plumeLoc.first, plumeLoc.second, 9000.0f); // 0, -12000
    
//    for(int s = 0; s < numskirts; s++)
//    {
//        skirt * sk = new skirt();
        
//        sk->init(center, 5000.0f, 200.0f, 15000.0f, 12000.0f, 100.0f, y);
//        // JG PARAM radius = 8000, step between samples = 200, target height = 15000, condensation height = 12000, thickness = 100
//        skirts.push_back(sk);
//        center.z() += 200.0f; // hardly any seperation
//    }

//    double previousStep = 0, d = 0;

//    PWM::Utils::writeTempImage("../ImageOutput/MergeTest/terrainModelTemp.png", y->getTerrain()->getTemperature());

//    std::cout << "Starting simulation..." << std::endl;

//    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
//    unsigned int plume_export_frequency = quickStepDT / plumeDT;
//    b.outputDirectory = "../ImageOutput/MergeTest/CapWilson";
//    for (int i = 0; i <= steps; ++i){
//        t.start();
        
//        plume.frame_draw(plumeDT);
      
//        for(int s = 0; s < numskirts; s++)
//        {
//            skirts[s]->uplift(plume, plumeDT, 8000.0f, 2000.0f); // JG PARAM - very broad influence zone horizontally, no tail off
//            skirts[s]->condensing();
//        }
     
//        //        transferDataFromPlumeToAtm(plume, b);
//        int frame_num = plume.frame_count / plume_export_frequency;
        
//        if (i % plume_export_frequency == 0)
//        {
//            plume.export_spheres_binary(plume_export_frequency, true);
//            for(int s = 0; s < numskirts; s++)
//                skirts[s]->exportSkirt(frame_num, s);
//        }

//        if (i % ((int) (quickStepDT / plumeDT)) == 0){
//            PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/CapWilson", i, cloudThreshold);
//            vdbExporter.write(b.getWorldModel(), "../output", frame_num, b.getSimTimePassed(), cloudThreshold); // uncomment for final run
//            PWM::Utils::writeAirImagesSelect("../ImageOutput/MergeTest/CapWilson", i, b.getWorldModel()); // images of atmo layers
//            if (i % ((int) (weatherDT / plumeDT)) != 0)
//                b.quick_internal_step();
//        }
//        if (i % ((int) (weatherDT / plumeDT)) == 0){
//            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
//            transferAshFromPlumeToAtm(plume, b.getWorldModel());

//            previousStep = b.getRunTimePassed();
//            b.step();
//            stripEdgeAsh(b.getWorldModel()->getAirLayers(), b.getWorldModel()->getConvectionLayers());
//            double currentStep = b.getRunTimePassed();
//            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
//            transferDataFromAtmToPlume(plume, b);
//        }
//        t.stop();
//        std::cerr << "Step " << i << " took " << t.peek() << "s" << std::endl;
//        if (i == steps / 8)
//            std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
//        else if (i == steps / 4)
//            std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
//        else if (i == (3 * steps / 8))
//            std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
//        else if (i == steps / 2)
//            std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
//        else if (i == (5 * steps / 8))
//            std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
//        else if (i == (3 * steps / 4))
//            std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
//        else if (i == (7 * steps / 8))
//            std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
//    }

//    b.printTiming(true);
//    for(int s = 0; s < numskirts; s++)
//        delete skirts[s];
//    return 0;
//}*/

// Ash rain variant that saves no data but performs looks for ashfall and cloud formation
int ashRainAnalysis(int steps, int width, int height, valType xSize, valType ySize, int layers, valType topOfAtmo){
    
    auto y = std::make_shared<PWM::Model::world<dsType, dsSType, valType, valType2>>(PWM::Model::world<dsType, dsSType, valType, valType2>());
    auto ter = constructTestTerrStruct("../resources/sthelens_detailed_sub5.obj", 400.f, xSize / 2);
    y->maxLayers = layers;
    y->atmoTop = topOfAtmo;    //Set up given heights for the air layers in the weather system
    std::vector<valType> layerHeights = {500, 1500, 2500, 2700, 3000, 3500, 4500, 5500, 6500};

    y->init("../resources/PlanetEarth.json", ter, width, height, xSize, ySize, layerHeights);

    const siv::PerlinNoise::seed_type tempSeed = 100u;
    const siv::PerlinNoise temPerNoiGen{ tempSeed };

    const siv::PerlinNoise::seed_type moisSeed = 647u;
    const siv::PerlinNoise moisPerNoiGen{ moisSeed };

    const siv::PerlinNoise::seed_type cloudSeed = 11883u;
    const siv::PerlinNoise cloudPerNoiGen{ cloudSeed };

    int octaves = 8;

    for (int k = 0; k < y->getAirLayers().size(); ++k){
        if (y->getAirLayer(k)->getHeight() < 7000){
            y->getAirLayer(k)->getVelocityTheta() = 1. + (0.005 * y->getAirLayer(k)->getHeight());
            y->getAirLayer(k)->getVelocityPhi() = 0;
        }
        else{
            y->getAirLayer(k)->getVelocityTheta() = 0;
            y->getAirLayer(k)->getVelocityPhi() = 1. + (0.005 * y->getAirLayer(k)->getHeight());
        }
        y->getAirLayer(k)->getCondensedWater() = 0.05;
        y->getAirLayer(k)->getMoisture() = 0.01;
        auto moisVal = (y->getAirLayer(k)->getHeight() < 3000) || (y->getAirLayer(k)->getHeight() > 8000) ? 0.05 : 0.12;
        auto cloudVal = (y->getAirLayer(k)->getHeight() < 4000) || (y->getAirLayer(k)->getHeight() > 10000) ? 0.02 : 0.08;
        for (int i = 0; i < height; ++i){
            for (int j = 0; j < width; ++j){
                valType newTemp = y->getAirLayer(k)->getTemperature(i, j) + (temPerNoiGen.normalizedOctave3D(i * 0.01, j * 0.01, k * 0.1, octaves) * 25);
                y->getAirLayer(k)->setTemperature(i, j, newTemp);
                valType newMois = moisPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * moisVal;
                y->getAirLayer(k)->setMoisture(i, j, newMois);
                valType newCloud = cloudPerNoiGen.normalizedOctave3D_01(i * 0.01, j * 0.01, k * 0.1, octaves) * cloudVal;
                y->getAirLayer(k)->setCondensedWater(i, j, newCloud);
            }
        }
    }

    float weatherDT = 2.f;
    float cveDT = 6.f;
    double plumeDT = 0.02;
    auto b = worldEngine(y, weatherDT, cveDT, true);
    
    scene_model plume;
    plume.setup_data(ter, y->getPlanet(), b.getWorldModel()->getAirLayers().size());
    plume.addAtmoModel(y);
    float ventHeight = ter.sampleHeightAt(0, 0);
    plume.setPlumeLoc(Eigen::Vector3f(0, 0, ventHeight));
    // take altitude layers from atmospheric model
    for (unsigned int i = 0; i<b.getWorldModel()->getAirLayers().size(); i++)
    {
        plume.wind_altitudes[i] = b.getWorldModel()->getAirLayer(i)->getHeight();
    }
    auto plumeLoc = convertFromCoordinatesToGridIdx(0, 0, y);
    std::cout << "Plume z_0 is at " << plume.z_0 << "m ASL and at grid location (" << plumeLoc.first << ", " << plumeLoc.second << ")." << std::endl;

    double previousStep = 0, d = 0;

    // Note : change plume_export_frequency so that the frame count export matches the weather export -> the files will be in /output
    unsigned int plume_export_frequency = 100;
    for (int i = 0; i <= steps; ++i){
        plume.frame_draw(plumeDT);
        int frame_num = plume.frame_count / plume_export_frequency;
        if (i % ((int) (weatherDT / plumeDT)) == 0){
            float currD = plume.run_time;
            transferDataFromPlumeToAtmAlt(plume, b.getWorldModel());
            transferAshFromPlumeToAtm(plume, b.getWorldModel());
            if (b.getStepCount() % 1 == 0){
                PWM::Utils::basicCloudExtract(b.getWorldModel(), "../ImageOutput/MergeTest/AshRain", i, 0.1);
//                PWM::Utils::airImagesAnalysis(i, 125, 125, 30, b.getWorldModel());
            }
            previousStep = b.getRunTimePassed();
            b.step();
            stripEdgeAsh(b.getWorldModel()->getAirLayers(), b.getWorldModel()->getConvectionLayers());
            transferDataFromAtmToPlume(plume, b);
            double currentStep = b.getRunTimePassed();
            std::cout << "Step " << b.getStepCount() << " of the weather engine took " << currentStep - previousStep << " seconds." << std::endl;
            std::cout << "Step " << b.getStepCount() << " of the plume took " << currD - d << " seconds." << std::endl;
            d = plume.run_time;
            if (i == steps / 8)
                std::cerr << "\033[1;36m12.5% of the way\033[0m" << std::endl;
            else if (i == steps / 4)
                std::cerr << "\033[1;36m25% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 8))
                std::cerr << "\033[1;36m37.5% of the way\033[0m" << std::endl;
            else if (i == steps / 2)
                std::cerr << "\033[1;36m50% of the way\033[0m" << std::endl;
            else if (i == (5 * steps / 8))
                std::cerr << "\033[1;36m62.5% of the way\033[0m" << std::endl;
            else if (i == (3 * steps / 4))
                std::cerr << "\033[1;36m75% of the way\033[0m" << std::endl;
            else if (i == (7 * steps / 8))
                std::cerr << "\033[1;36m87.5% of the way\033[0m" << std::endl;
        }
    }
    int slCount = plume.smoke_layers.size() - 1;
    std::cout << "One end of the plume is at " << plume.smoke_layers[0].center.z() << "m and moving at v = " << plume.smoke_layers[0].v << " m.s-1" << std::endl;
    std::cout << "The other end of the plume is at " << plume.smoke_layers[slCount].center.z() << "m and moving at v = " << plume.smoke_layers[slCount].v << " m.s-1" << std::endl;

    b.printTiming(true);
    return 0;
}

int main(int argc, char** argv){
    if (argc > 1){
        if (strcasecmp(argv[1], "DC") == 0 || strcasecmp(argv[1], "Donut") == 0 || strcasecmp(argv[1], "DonutClouds") == 0)
            return donutClouds(12500, 512, 512, 30000.0, 30000.0, 20, 20000);
        else if (strcasecmp(argv[1], "AR") == 0 || strcasecmp(argv[1], "Ash") == 0 || strcasecmp(argv[1], "AshRain") == 0)
            return ashRain(12500, 512, 512, 30000.0, 30000.0, 20, 20000);
        else if (strcasecmp(argv[1], "ARJG") == 0 || strcasecmp(argv[1], "AshJG") == 0 || strcasecmp(argv[1], "AshRainJG") == 0)
            return ashRain(12500, 512, 512, 30000.0, 30000.0, 20, 20000, true);
        else if (strcasecmp(argv[1], "DCBG") == 0 || strcasecmp(argv[1], "DonutBig") == 0 || strcasecmp(argv[1], "DonutCloudsBigGrid") == 0)
            return donutCloudsBigScale(12500, 20, 20000);
        else if (strcasecmp(argv[1], "ARBG") == 0 || strcasecmp(argv[1], "AshBig") == 0 || strcasecmp(argv[1], "AshRainBigGrid") == 0)
            return ashRainBigScale(180000, 20, 20000);
        else if (strcasecmp(argv[1], "ARBGJG") == 0 || strcasecmp(argv[1], "AshBigJG") == 0 || strcasecmp(argv[1], "AshRainBigGridJG") == 0)
            return ashRainBigScale(180000, 20, 20000, true, false, true);
        else if (strcasecmp(argv[1], "ARBGA") == 0 || strcasecmp(argv[1], "AshBigAnalysis") == 0 || strcasecmp(argv[1], "AshRainBigGridAnalysis") == 0)
            return ashRainBigScale(180000, 20, 20000, false, true);
        else if (strcasecmp(argv[1], "ARBGF") == 0 || strcasecmp(argv[1], "AshBigFinal") == 0 || strcasecmp(argv[1], "AshRainBigGridFinal") == 0)
            return ashRainBigScaleFinal(60000, 20, 20000, 100); // 300
        else if (strcasecmp(argv[1], "ARV") == 0 || strcasecmp(argv[1], "ARVAL") == 0 || strcasecmp(argv[1], "AshRainVal") == 0 || strcasecmp(argv[1], "AshRainValidation") == 0)
            return ashRainValidation(180000, 20, 20000, 512, 75000, 300);
        else if (strcasecmp(argv[1], "S") == 0 || strcasecmp(argv[1], "Skirt") == 0)
            return skirts3(32000, 20, 20000);
        else if (strcasecmp(argv[1], "CW") == 0 || strcasecmp(argv[1], "CapWilson") == 0)
            return capWilson(6000, 20, 12000);
        else{
            std::cerr << "Incorrect option chosen! Options are DC, AR, ARJG, DCBG, ARBG, ARBGJG, ARBGA, ARBGF, S, or CW!" << std::endl;
            return 1;
        }
    }
    else
        return capWilson(32000, 20, 20000);
}
