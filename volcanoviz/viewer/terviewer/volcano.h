/*******************************************************************************
 *
 * EcoSynth - Data-driven Authoring of Large-Scale Ecosystems (Undergrowth simulator)
 * Copyright (C) 2020  J.E. Gain  (jgain@cs.uct.ac.za)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 ********************************************************************************/


/**
 * @file
 */

#ifndef VOLCANO_H
#define VOLCANO_H

#include <memory>
#include "vecpnt.h"
#include "view.h"
#include "trenderer.h"
#include "common/basic_types.h"
#include "common/cnpy.h"

// volcano.h: for managing volcano and weather data
// author: James Gain
// date: 17 December 2012

enum class VolcanoMode
{
    TEMPERATURE,
    MOISTURE,
    CONVECTION,
    ASH,
    VELOCITYX,
    VELOCITYY,
    CONDENSED,
    VMEND
};

const std::array<VolcanoMode, 7> all_volcanomodes = {
    VolcanoMode::TEMPERATURE,
    VolcanoMode::MOISTURE,
    VolcanoMode::CONVECTION,
    VolcanoMode::ASH,
    VolcanoMode::VELOCITYX,
    VolcanoMode::VELOCITYY,
    VolcanoMode::CONDENSED
}; // to allow iteration over the brushtypes

/// Storing per-frame collections of spheres


struct PlumeSphere
{
    vpPoint pos;
    float radius;
    float density;
    float charge;
    int sphtype;
};

struct AirLayer
{
    float altitude;
    bool active;
};

class AirLayers
{

private:
    std::vector<AirLayer> planes;
    std::vector<std::vector<QImage>> temperature;
    std::vector<std::vector<QImage>> moisture;
    std::vector<std::vector<QImage>> convection;
    std::vector<std::vector<QImage>> ash;
    std::vector<std::vector<QImage>> velocityx;
    std::vector<std::vector<QImage>> velocityy;
    std::vector<std::vector<QImage>> condensed;
    int numframes, currframe, numlayers, currlayer;
    VolcanoMode vmode;
    
public:

    /// Constructor
    AirLayers()
    {
        numframes = 0;
        currframe = 0;
        numlayers = 0;
        currlayer = 0;
        vmode = VolcanoMode::TEMPERATURE;
    }

    /// Destructor
    ~AirLayers()
    {
        clear();
      
    }
    
    /**
     * clear: empty the plume sphere list
     */
    void clear()
    {
        planes.clear();
        for(auto &t: temperature)
            t.clear();
        for(auto &m: moisture)
            m.clear();
    }
    
    // getter and setter for currently active frame, layer, mode, etc.
    void setNumFrames(int numf){ numframes = numf; }
    int getNumFrames(){ return numframes; }
    int getNumLayers(){ return numlayers; }
    
    int getFrame(){ return currframe; }
    void setFrame(int f)
    {
        if(f >= 0 && f < numframes)
            currframe = f;
        else
            std::cerr << "Error AirLayers::setFrame - frame " << f << " out of range" << std::endl;
    }
    
    int getLayer(){ return currlayer; }
    void setLayer(int lyr)
    {
        if(lyr >= 0 && lyr < numlayers)
            currlayer = lyr;
        else
            std::cerr << "Error AirLayers::setLayer - layer " << lyr << " out of range" << std::endl;
    }
    
    void setMode(VolcanoMode vizmode){ vmode = vizmode; }
    void setActiveStatus(bool status){ planes[currlayer].active = status; }
    bool getActiveStatus(){ return planes[currlayer].active; }
    
    std::vector<AirLayer> & getCurrentPlanes(){ return planes; }
    QImage & getCurrentLayerImage();
    QImage & getCurrentMoisture(){ return moisture[currframe][currlayer]; }
    QImage & getCurrentTemperature(){ return temperature[currframe][currlayer]; }
    QImage & getCurrentConvection(){ return convection[currframe][currlayer]; }
    QImage & getCurrentAsh(){ return ash[currframe][currlayer]; }
    QImage & getCurrentVelocityX(){ return velocityx[currframe][currlayer]; }
    QImage & getCurrentVelocityY(){ return velocityy[currframe][currlayer]; }
    QImage & getCurrentCondensed(){ return condensed[currframe][currlayer]; }

    /**
       * Load airlayers from file.
       * @param filename   File to load (simple ascii format)
       * @param interval   Number of timesteps between airlayer writes
       */
    void loadAirLayers(string &basefile, int interval = 1);
};

/// Sample air parcel in the skirt / cap
struct Parcel
{
    vpPoint position;         ///< 3d spatial position
    float vapour;             ///< water vapour value in kg / m^3
    float condensation;       ///< condensed water value in kg / m^3
};

struct Skirt
{
    std::vector<std::vector<Parcel>> frames; ///< per frame parcels making up the deformed air layer
    float thickness;                    ///< thickness of the skirt in meters
};

/// a single thin cloud layer subject to plume uplift
class Skirts
{
private:
  
    std::vector<Skirt> skirtlayers; ///< per layer skirts
    int numframes, currframe;       ///< frame control
    int numlayers;                  ///< number of stacked skirts
    
public:
    
    Skirts()
    {
        numframes = 0; currframe = 0; numlayers = 0;
    }
    
    ~Skirts()
    {
        clear();
    }
    
    // empty the layer data
    void clear()
    {
        for(auto &s: skirtlayers)
            for(auto &f: s.frames)
                f.clear();
    }
    
    // getter and setter for currently active frame
    void setFrame(int f)
    {
        if(f >= 0 && f < numframes)
            currframe = f;
        else
            std::cerr << "Error skirts::setFrame - frame " << f << " out of range" << std::endl;
    }
    
    void setNumFrames(int numf){ numframes = numf; }
    int getNumFrames(){ return numframes; }
    int getFrame(){ return currframe; }
    int getNumLayers(){ return numlayers; }
    
    /// Obtain grid size @a dx and @a dy
    // void getGridDim(int & dx, int & dy) const;
    int getNumSamples(int l);
    
    /// get the position and condensation value for a particular sample
    void getFrameSample(int l, int s, vpPoint &pos, float &cval);
        
    /// read sequence of skirts from files
    void loadSkirts(string &basefile, int numskirts, int interval);
};

class Plume
{

private:
    std::vector<std::vector<PlumeSphere>> plumespheres;
    int numframes, currframe;
    
public:

    /// Constructor
    Plume()
    {
        numframes = 0;
        currframe = 0;
    }

    /// Destructor
    ~Plume()
    {
        clear();
      
    }
    
    /**
     * clear: empty the plume sphere list
     */
    void clear()
    {
        for(auto &f: plumespheres)
            f.clear();
    }
    
    // getter and setter for currently active frame
    void setFrame(int f)
    {
        if(f >= 0 && f < numframes)
            currframe = f;
        else
            std::cerr << "Error Plume::setFrame - frame " << f << " out of range" << std::endl;
    }
    
    void setNumFrames(int numf){ numframes = numf; }
    int getNumFrames(){ return numframes; }
    
    int getFrame(){ return currframe; }
    
    std::vector<PlumeSphere> & getCurrentPlume(){ return plumespheres[currframe]; }

    // draw: display terrain using OpenGL as a triangle mesh
    void draw(View * view, PMrender::TRenderer *renderer) const;

    /**
       * Load plume spheres from file.
       * @param filename   File to load (simple ascii format)
       */
    void loadPlume(string &basefile, int start);
};

#endif // VOLCANO_H
