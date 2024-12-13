#include "skirt.hpp"
#include "mathUtils.h"

#include <sstream>
#include <streambuf>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

/*
void skirt::getGridDim(int & dx, int & dy) const
{
    dx = height->width();
    dy = height->height();
}

void skirt::convertGridToWorld(int gx, int gy, float &wx, float &wy, float &wz)
{
    // find position in world coordinate space
    int dx, dy;
    float convx, convy, h;

    getGridDim(dx, dy);
    convx = extent.x() / (float) (dx-1);
    convy = extent.y() / (float) (dy-1);
    h = height->get(gx, gy);
    
    wx = start.x() + (float) gx * convx;
    wy = start.y() + (float) gy * convy;
    wz = start.z() + h;
}

void skirt::transferWaterFromAtmosphere(int gx, int gy)
{
    float cval, vval;
    float x, y, z;
    
    // find closest water vapour and condensation values in atmospheric layer
    // consider adding small level of noise
 
    // convert location to atmospheric coordinates
    convertGridToWorld(gx, gy, x, y, z);
    
    // lookup values
    vval = atmoModel->getLocalMoist(x, y, z);
    cval = atmoModel->getLocalClouds(x, y, z);
    
    condensation->set(gx, gy, cval);
    vapour->set(gx, gy, vval);
 }
 */
   
void skirt::uplift(scene_model &plume, float dt, float hzone, float vzone)
{
    float lift;
    float maxhght = 0.0f;
    int substeps;
    auto lapseRate = this->atmoModel->getPlanet()->getTempLapseRate() * -1;
    
    if(!targetreached)
    {
        plume.boundSpheres(vzone);
        substeps = plume.calcSubsteps();
        #pragma omp parallel for
        for(auto &s: samples)
        {
            lift = plume.calcSkirtUplift(s.position, dt, hzone, vzone, substeps);
            s.position.z() += lift;
            s.temp += lift * lapseRate;
        }
        #pragma omp barrier

        #pragma omp parallel for reduction(max:maxhght)
        for(auto &s: samples)
        {
            maxhght = (s.position.z() > maxhght) ? s.position.z() : maxhght;
        }
        #pragma omp barrier

        targetreached = (maxhght >= targetheight);
        if(targetreached)
            std::cerr << "SKIRT TARGET REACHED" << std::endl;
    }
    
        
/*
    for(int gx = 0; gx < dx; gx++)
        for(int gy = 0; gy < dy; gy++)
        {
            convertGridToWorld(gx, gy, wx, wy, wz);
            gridpoint = Eigen::Vector3f(wx, wy, wz);
            
            lift = plume.calcSkirtUplift(gridpoint);
            height->set(gx, gy, height->get(gx, gy) + lift);
            if(lift > maxlift)
                maxlift = lift;
            avglift += lift;
        }
 */
    // cerr << "Maximum lift = " << maxlift << " average lift = " << avglift / ((float) dx*dy) << endl;
}

void skirt::phaseTransition()
{
    float wx, wy, wz, maxwater = 0.0f;
    Eigen::Vector3f gridpoint;

    #pragma omp parallel for
    for(auto &s: samples)
    {
        wx = s.position.x(); wy = s.position.y(); wz = s.position.z();
        
        // values required for condensation calculation
        float temp = s.temp;
        float pres = s.pressure;
        float wilsonMod = 0;
        if (pres != 0){
            if (pres < 0)
                wilsonMod = std::abs(pres) * 0.5;
            else{
                temp *= 1 / (1 + pres);
                wilsonMod = pres;
            }
        }
        float moisture = s.vapour;
        float satVapPres = PWM::Utils::saturationVapourPressure(temp);
        float altitude = wz;
        float atmPres = PWM::Utils::altitudeAdjustedPressure(altitude, atmoModel->getPlanet());
        float density = PWM::Utils::pressureAdjustedDensity(atmPres, temp, atmoModel->getPlanet());
        float humidity = moisture / density;
        float partialPressure = atmPres * (humidity / (0.622 + humidity));
        float satHumidity = (0.622 * satVapPres) / (atmPres - satVapPres);

        // use standard condensation procedure, exluding heat exchange
           
        if (moisture > 0.0f && (satVapPres < partialPressure) && (atmPres > satVapPres))
        {
            float diff = humidity - satHumidity;
            float water = diff * 0.1f  * (1 + wilsonMod); // condensationCoefficient;

            if (moisture < water)
                water = moisture;
            s.vapour -= water;
            s.condensation += water;
            if(water > maxwater)
                maxwater = water;
        }
        else if ((s.condensation > 0) && (satVapPres > partialPressure) && (atmPres > satVapPres)){
            float diff = satHumidity - humidity;
            float water = diff * 0.1f; // evaporationCoefficient
            if (s.condensation < water)
                water = s.condensation;
            s.condensation -= water;
            s.vapour += water;
        }
    }
    #pragma omp barrier
    // cerr << "Max skirt condensation = " << maxwater << endl;
    
}

void skirt::updateTemp()
{
#pragma omp parallel for
    for(auto &s: samples)
    {
        float wx = s.position.x(), wy = s.position.y(), wz = s.position.z();
        // values required for condensation calculation
        float temp = atmoModel->getLocalTemp(wx, wy, wz);
    }
#pragma omp barrier
}

void skirt::init(Eigen::Vector3f i_center, float i_rad, float i_step, float i_targethght, float i_condenseThres, float i_thickness, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& i_model)
{
    center = i_center; radius = i_rad; step = i_step; targetheight = i_targethght; condensationThres = i_condenseThres; thickness = i_thickness;
    atmoModel = i_model;
    
    // establish number of gridpoints for square domain based on the step size
    int dx, dy;
    auto gridLen = i_model->getTerrain()->getElevation().gridLength();
    
    Eigen::Vector3f start = Eigen::Vector3f(center.x() - radius, center.y() - radius, center.z());
    float extent = 2.0f * radius;
    float steps = (int) (extent / step);
    targetreached = (center.z() >= targetheight);
   
    for(int x = 0; x < steps; x++)
        for(int y = 0; y < steps; y++)
        {
            parcel p;
            p.position = Eigen::Vector3f(start.x() + x * i_step, start.y() + y * i_step, start.z());
            float i = (i_model->getWorldCellWidth() / 2) + x/gridLen;
            float j = (i_model->getWorldCellWidth() / 2) + y/gridLen;
            double terHeight = i_model->getTerrain()->getElevation().sampleAt(i, j);
            if (p.position.z() <= terHeight)
                continue;
            
            if((p.position - center).norm() < radius)
            {
                // transfer water from atmosphere
                // alternative is to seed it randomly
                p.vapour = atmoModel->getLocalMoist(p.position.x(), p.position.y(), p.position.z());
                p.condensation = atmoModel->getLocalClouds(p.position.x(), p.position.y(), p.position.z());
                p.temp = atmoModel->getLocalTemp(p.position.x(), p.position.y(), p.position.z());
                p.pressure = 0;
                p.timeSinceExpl = 0;
                samples.push_back(p);
            }
        }
}

void skirt::exportSkirt(int framenum, int layernum)
{
    string skirtfname, framenumstr;
   
    framenumstr = to_string(framenum);
    int precision = 3 - std::min((int) framenumstr.size(), 3);
    string framename = std::string(precision, '0').append(framenumstr);
    skirtfname = "../output/skirts_" + framename + "_" + to_string(layernum) + ".txt";
    write(skirtfname);
}

void skirt::write(const std::string &filename)
{
    ofstream outfile;

    outfile.open((char *) filename.c_str(), ios_base::out);
    if(outfile.is_open())
    {
        // write header
        outfile << (int) samples.size() << " " << thickness / 1000.0f << endl;
       
        // write out samples
        for(auto & s: samples)
        {
            if(s.condensation > condensationThres) // does it meet density threshold to be visible?
                outfile << s.position.x()/1000.0f << " " << s.position.z()/1000.0f << " " << s.position.y()/1000.0f << " " << s.condensation << endl;
            
        }
        outfile.close();
    }
    else
    {
        cerr << "Error skirt::write - unable to open file " << filename << endl;
    }   
}

skirt::skirt(Eigen::Vector3f i_center, float i_rad, float i_step, float i_targethght, float i_condenseThres, float i_thickness, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2> > &i_model)
{
    this->init(i_center, i_rad, i_step, i_targethght, i_condenseThres, i_thickness, i_model);
}

wilsonCloud::wilsonCloud(Eigen::Vector3f i_center, float i_radius, float i_step, float startTime, float lifeTime, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2> > &i_model) : skirt(i_center, i_radius, i_step, i_center.z(), 0.1, i_radius, i_model)
{
    this->samples = std::vector<parcel>();
    expStart = startTime;
    expLife = lifeTime;
    this->init(i_center, i_radius, i_step, i_center.z() + 800, 0.1, i_radius, i_model);
}

wilsonCloud::~wilsonCloud()
{
    samples.clear();
}

void wilsonCloud::updateExp(float dT)
{
    if (explosion == nullptr)
        return;

    explosion->updateWave(dT);

    //If explosion isn't started yet, do nothing.
    if (!explosion->isStarted()){
        return;
    }
    this->expLife += dT;

    float expMaxLife = explosion->getMaxLife();
    reset = this->expLife > (expMaxLife + 10 * dT);

    float zero = 0;
    if (reset){
        resetOnce = true;
        return;
    }

    if (!resetOnce){
        //Calculate limited area for optimisation
        auto worldLoc = explosion->getCentre();
        float expCentreHeight = explosion->getHeight();
        float radius = explosion->getRadius();
        float expPower = explosion->getInitialPower();
        auto eLoc = std::make_tuple(worldLoc.first, worldLoc.second, expCentreHeight);

        #pragma omp parallel for
        for (auto i : this->samples){
            auto loc = std::make_tuple(i.position.x(), i.position.y(), i.position.z());
            float dist = PWM::Utils::calcCartesianDistance<float>(eLoc, loc);
            if (dist < radius){
                float t = i.timeSinceExpl;
                float pres = i.pressure;
                if (pres != 0){
                    i.timeSinceExpl = t + dT;
                    float initVal = pres * (t + 1);
                    float newVal = initVal / std::pow(i.timeSinceExpl + 1, 1);
                    if (t > pressureWaveLife)
                        newVal *= -1;
                    i.pressure = newVal;
                }
                else{
                    i.pressure = expPower / std::pow(dist, 1);
                }
            }
        }
        #pragma omp barrier
    }
}

void wilsonCloud::addExplosion(std::shared_ptr<PWM::Model::pressureWave<double> > &exp)
{
    this->explosion = exp;
}

void wilsonCloud::uplift(scene_model &plume, float dt, float hzone, float vzone)
{
    return;
}

void wilsonCloud::init(Eigen::Vector3f i_center, float i_rad, float i_step, float i_targethght, float i_condenseThres, float i_thickness, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2> > &i_model)
{
    // establish number of gridpoints for square domain based on the step size
    int dx, dy;

    Eigen::Vector3f start = Eigen::Vector3f(i_center.x() - i_rad, i_center.y() - i_rad, i_center.z() - i_rad);
    std::cerr << "Start position of (" << start.x() << ", " << start.y() << ", " << start.z() << ")" << std::endl;
    float extent = 2.0f * i_rad;
    int steps = std::ceil(extent / i_step);

    samples.clear();
    auto gridLen = i_model->getTerrain()->getElevation().gridLength();

    for(int x = 0; x < steps; x++)
        for (int z = 0; z < steps; z++){
            for(int y = 0; y < steps; y++)
            {
                parcel p;
                p.position = Eigen::Vector3f(start.x() + x * i_step, start.y() + y * i_step, start.z() + z * i_step);
                float i = (i_model->getWorldCellWidth() / 2) + x/gridLen;
                float j = (i_model->getWorldCellWidth() / 2) + y/gridLen;
                double terHeight = i_model->getTerrain()->getElevation().sampleAt(i, j);
                if (p.position.z() <= terHeight)
                    continue;

                if((p.position - i_center).norm() < i_rad)
                {
                    // transfer water from atmosphere
                    // alternative is to seed it randomly
                    p.vapour = i_model->getLocalMoist(p.position.x(), p.position.y(), p.position.z());
                    p.condensation = i_model->getLocalClouds(p.position.x(), p.position.y(), p.position.z());
                    p.temp = i_model->getLocalTemp(p.position.x(), p.position.y(), p.position.z());
                    p.pressure = 0;
                    p.timeSinceExpl = 0;
                    samples.push_back(p);
                }
            }
    }
}
