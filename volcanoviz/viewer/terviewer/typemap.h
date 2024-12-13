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


// TypeMap class for passing around terrian type information

#ifndef _TYPEMAP
#define _TYPEMAP

#include "glheaders.h"
#include <vector>
#include <memory>
#include <common/region.h>
#include "common/basic_types.h"
#include <array>

#define def_gentleslope 0.3f
#define def_steepslope 0.7f

const float gentlecol[] = {1.0f, 0.631f, 0.521f, 1.0f};
const float steepcol[] = {0.671f, 0.331f, 0.221f, 1.0f};
const float cliffcol[] = {0.271f, 0.031f, 0.0f, 1.0f};

// const float freecol[] = {0.874f, 0.827f, 0.683f, 1.0f};
const float freecol[] = {0.7529f, 0.7137f, 0.6549f, 1.0f};
const float watercol[] = {0.1098f, 0.522f, 0.741f, 1.0f};
// const float icecol[] = {0.812f, 0.933f, 0.98f, 1.0f};
const float icecol[] = {0.753f, 0.87f, 0.98f, 1.0f};
const float rockcol[] = {0.447f, 0.450f, 0.496f, 1.0f};
const float trailcol[] = {0.894f, 0.769f, 0.531f, 1.0f};
const float grasscol[] = {0.573f, 0.757f, 0.565f, 1.0f};
const float shrubscol[] = {0.310f, 0.416f, 0.310f, 1.0f};
const float forestcol[] = {0.289f, 0.530f, 0.407f, 1.0f};
const float snowcol[] = {1.0f, 1.0f, 1.0f, 1.0f};
const float obstaclecol[] = {0.850f, 0.465f, 0.461f, 1.0f};

// const float boids1col[] = {0.545f, 0.741f, 0.537f, 1.0f};
const float boids1col[] = {0.682f, 0.902f, 0.474f, 1.0f}; // 6
const float boids2col[] = {0.756f, 1.0f, 0.315f, 1.0f}; // 8
const float boids3col[] = {0.573f, 0.757f, 0.258f, 1.0f}; // 12
const float boids4col[] = {0.447f, 0.588f, 0.439f, 1.0f}; // 1
const float boids5col[] = {0.349f, 0.463f, 0.345f, 1.0f}; // 3
const float boids6col[] = {0.400f, 0.757f, 0.314f, 1.0f}; // 14
const float boids7col[] = {0.390f, 0.653f, 0.379f, 1.0f};
const float boids8col[] = {0.374f, 0.644f, 0.362f, 1.0f};
// const float boids9col[] = {1.0f, 0.71f, 0.76f, 1.0f};
// const float boids10col[] = {1.0f, 0.71f, 0.76f, 1.0f};

const float boidsE1col[] = {0.801f, 0.827f, 0.682f, 1.0f};
const float boidsE2col[] = {0.875f, 0.772f, 0.682f, 1.0f};
const float boidsE3col[] = {0.925f, 0.769f, 0.533f, 1.0f};
const float boidsE4col[] = {0.774f, 0.754f, 0.600f, 1.0f};
const float boidsE5col[] = {0.808f, 0.788f, 0.651f, 1.0f};
const float boidsE6col[] = {0.910f, 0.816f, 0.510f, 1.0f};
const float boidsP1col[] = {0.801f, 0.827f, 0.682f, 1.0f};
const float boidsP2col[] = {0.875f, 0.772f, 0.682f, 1.0f};
const float boidsP3col[] = {0.925f, 0.769f, 0.533f, 1.0f};
const float boidsP4col[] = {0.774f, 0.754f, 0.600f, 1.0f};
const float boidsP5col[] = {0.808f, 0.788f, 0.651f, 1.0f};
const float boidsP6col[] = {0.910f, 0.816f, 0.510f, 1.0f};
const float addcol[] = {0.75f, 0.75f, 0.75f, 1.0f};
const float delcol[] = {0.25, 0.25f, 0.25f, 1.0f};

enum class BrushType
{
    FREE,                  //< neutral type
    ICE,                   //< hard frozen surface
    WATER,                 //< lakes, rivers, and sea
    ROCK,                  //< no plant life
    TRAIL,                 //< preferential worn area where herds have passed repeatedly
    GRASS,                 //< freely passable area with possiblity for grazing
    SHRUBS,                //< difficult but not impossible to pass through, possible grazing
    FOREST,                //< denser, less passable area
    SNOW,                  //< deep soft snow, slows progress
    OBSTACLE,              //< impassable area such as a fence
    // note that the last 3 paint to a seperate slope map
    GENTLE,                //< no restrictions on passage
    STEEP,                 //< slope is above a threshold that slows progress
    CLIFF,                 //< impassable terrain for most species
    BOIDS1,                //< slots for boid behaviour brushes
    BOIDS2,
    BOIDS3,
    BOIDS4,
    BOIDS5,
    BOIDS6,
    HERDADD,               //< randomly place herd animals in brush region
    HERDDEL,               //< delete herd animals from brush region
    BLTEND
};

enum class TypeMapType
{
    DISPLAY,        //< map for user display, combining other maps
    BOIDS,          //< to show boid behaviour categories
    CATEGORY,       //< to show terrain attribute categories
    SLOPE,          //< to show slope categories
    TRAIL,          //< to show an enclosing shape around herd movement
    HOME,           //< starting area for boids
    TMTEND
};
const std::array<TypeMapType, 6> all_typemaps = {TypeMapType::DISPLAY, TypeMapType::BOIDS, TypeMapType::CATEGORY, TypeMapType::SLOPE, TypeMapType::TRAIL, TypeMapType::HOME}; // to allow iteration over the typemaps

class TypeMap
{
private:
    MapInt * tmap;      ///< a map corresponding to the terrain storing integer types
    std::vector<GLfloat *> colmap;  ///< a 32-element lookup table for converting type indices to colours
    Region dirtyreg;                ///< bounding box in terrain grid integer coordinates (e.g, x=[0-width), y=[0-hieght))
    TypeMapType usage;              ///< indicates map purpose
    int numSamples;                 ///< number of active entries in lookup table

    /// Set up the colour table with colours appropriate to the initial ecosystem pallete of operations
    void initPaletteColTable();

    /**
     * @brief initPerceptualColTable Set up a colour table sampled from a perceptually uniform colour map stored in a CSV file
     * @param colmapfile        file on disk containing the colour map
     * @param samples           number of samples taken from the colour map
     * @param truncend          proportion of the colourmap to select, truncating from the upper end
     */
    void initPerceptualColTable(std::string colmapfile, int samples, float truncend = 1.0f);

    /// clip a region to the bounds of the map
    void clipRegion(Region &reg);
    
public:

    TypeMap(){ usage = TypeMapType::DISPLAY; }

    TypeMap(TypeMapType purpose, int fillval = -1);

    /**
     * Create type map that matches the terrain dimensions
     * @param w         map width
     * @param h         map height
     * @param purpose   map purpose, to represent different kinds of layers
     * @param fillval   category value applied to the entire map
     */
    TypeMap(int w, int h, TypeMapType purpose, int fillval = -1);

    virtual ~TypeMap();

    /// getters for width and height
    int width(){ return tmap->width(); }
    int height(){ return tmap->height(); }

    /// fill map with a certain colour
    void fill(int val){ tmap->fill(val); }

    /// Match type map dimensions to @a w (width) and @a h (height)
    void matchDim(int w, int h);
    
    /// clear typemap to unconstrained
    void clear();

    /// getter for underlying map
    MapInt * getMap(void){ return tmap; }
    
    /// getter for individual value
    int get(int x, int y){ return tmap->get(y,x); }
    void set(int x, int y, int val){ tmap->set(y,x,val); }
    
    /// replace underlying map
    void replaceMap(MapInt * newmap);

    /// load from file, return number of clusters
    int load(const std::string &filename, TypeMapType purpose);

    /// load category data from PNG file, return number of clusters
    bool loadCategoryImage(const std::string &filename);

    /// load category data from text file
    bool readCategoryFile(const std::string &filename);

    /// convert a floating point map into a discretized type map
    int convert(MapFloat * map, TypeMapType purpose, float range);

    /// save a mask file version of the type map
    void save(const std::string &filename);

    /**
     * @brief saveToPainrImage   Save paint map out as a greyscale image
     * @param filename      Name of file to save to
     */
    void saveToPaintImage(const std::string &filename);

    /// getter for colour table
    std::vector<GLfloat *> * getColourTable(void) { return &colmap; }

    /// getter for sampling range
    int getTopSample(){ return numSamples-1; }

    /// getter for update region
    Region getRegion(void) { return dirtyreg; }

    /// setter for update region
    void setRegion(const Region& toupdate)
    {
        dirtyreg = toupdate;
        clipRegion(dirtyreg);
    }

    /// return region that covers the entire type map
    Region coverRegion()
    {
        return Region(0,0,width(),height());
    }

    /// setter for purpose
    void setPurpose(TypeMapType purpose);

    /// Reset the indicated type to zero everywhere it appears in the map
    void resetType(int ind);

    /**
     * Index to colour translation. First element is the erase colour
     * @param ind   index for a type (potentially) stored in the type map
     * @retval      4-element RGBA colour associated with @a ind
     */
    GLfloat * getColour(int ind)
    {
        if(ind >= 0 && ind < 32) // limit on number of available colours
            return colmap[ind];
        else
            return NULL;
    }

    /**
     * Set the colour associated with a particular index
     * @param ind   index for the type that is to be given a new colour association
     * @param col   4-element RBGA colour to associate with @a ind
     */
    void setColour(int ind, GLfloat * col);

    int getNumSamples();

    template<typename T>
    int convert(const T &map, TypeMapType purpose, float range)
    {
        int tp, maxtp = 0;
        int width, height;

        map.getDim(width, height);
        matchDim(width, height);

        for(int x = 0; x < width; x++)
            for(int y = 0; y < height; y++)
            {
                tp = 0;
                switch(purpose)
                {
                    case TypeMapType::BOIDS: // do nothing
                        break;
                    case TypeMapType::CATEGORY: // do nothing, since categories are integers not floats
                        break;
                    /*
                    case TypeMapType::SLOPE:
                        val = map.get(x, y);
                        if(val > maxval)
                            maxval = val;

                        // discretise into ranges of illumination values
                        // clamp values to range
                        if(val < 0.0f) val = 0.0f;
                        if(val > range) val = range;
                        tp = (int) (val / (range+_pluszero) * (numSamples-2)) + 1;
                        break;
                    */
                    default:
                        break;
                }
                tmap->set(y,x,tp);

                if(tp > maxtp)
                    maxtp = tp;
            }
        /*
        if(purpose == TypeMapType::CDM)
        {
            cerr << "Minimum colour value = " << mincm << endl;
            cerr << "Maxiumum colour value = " << maxcm << endl;
        }*/
        return maxtp;
    }
};

#endif
