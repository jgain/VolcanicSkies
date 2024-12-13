/**
 description: Sampled planar layers uplifted by plume and used to represent bell and skirt clouds
 author: James Gain
 date: July 2023
 */

#ifndef CF_SKIRT_HPP
#define CF_SKIRT_HPP



#include "planet.h"
#include "smoke.hpp"
#include "smokeLayer.hpp"
#include "terrain_structure.hpp"
#include "world.h"
#include <string>

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;

/// grid of regularly sampled values
class MapFloat
{
private:
    int gx, gy;                     //< grid dimensions
    std::vector<float> fmap;        //< grid of floating point values

public:

    MapFloat(){ gx = 0; gy = 0; initMap(); }

    ~MapFloat(){ delMap(); }

    /// return the row-major linearized value of a grid position
    inline int flatten(int dx, int dy) const { return dy * gx + dx; }

    inline void idx_to_xy(int idx, int &x, int &y) const
    {
        x = idx % gx;
        y = idx / gx;
    }

    /// getter for grid dimensions
    void getDim(int &dx, int &dy) const
    { dx = gx; dy = gy; }

    int height() const { return gy; }
    int width() const{ return gx; }

    /// setter for grid dimensions
    void setDim(int dx, int dy){ gx = dx; gy = dy; initMap(); }
    template<typename T>
    void setDim(const T &other)
    {
        int w, h;
        other.getDim(w, h);
        this->setDim(w, h);
    }

    template<typename T>
    void clone(const T &other)
    {
        setDim(other);
        int w, h;
        other.getDim(w, h);
        for (int y = 0; y < h; y++)
        {
            for (int x = 0; x < w; x++)
            {
                set(x, y, other.get(x, y));
            }
        }
    }

    /// clear the contents of the grid to empty
    void initMap(){ fmap.clear(); fmap.resize(gx*gy); }

    /// completely delete map
    void delMap(){ fmap.clear(); }

    /// set grass heights to a uniform value
    void fill(float h){ fmap.clear(); fmap.resize(gx*gy, h); }

    /// getter and setter for map elements
    float get(int x, int y) const { return fmap[flatten(x, y)]; }
    float get(int idx) const{ return fmap[idx]; }
    void set(int x, int y, float val){ fmap[flatten(x, y)] = val; }

    /// get pointer to the raw map structure
    float * getPtr(){ return &fmap[0]; }

    /**
     * @brief read  read a floating point data grid from file
     * @param filename  name of file to be read
     * @return true if the file is found and has the correct format, false otherwise
     */
    bool read(std::string filename);

    float *data() { return fmap.data(); }

    std::vector<float>::iterator begin()
    {
        return fmap.begin();
    }

    std::vector<float>::iterator end()
    {
        return fmap.end();
    }

    /**
     * @brief getBilinear returns the bilinearly interpolated value between four cells
     * @param x  coordinate in range [0..gx-1]
     * @param y  coordinate in range [0..gy-1]
     * @return bilinearly interpolated value between floor(x,y) and ceil(x,y)
     */
    float getBilinear(float x, float y) const
    {
        int i0 = std::max(0, std::min(gx-1, int(x)));
        int j0 = std::max(0, std::min(gy-1, int(y)));
        int i1 = std::max(0, std::min(gx-1, i0+1));
        int j1 = std::max(0, std::min(gy-1, j0+1));
        float fx = std::max(0.0f, x - float(i0));
        float fy = std::max(0.0f, y - float(j0));

        float v00 = get(i0, j0);
        float v01 = get(i0, j1);
        float v10 = get(i1, j0);
        float v11 = get(i1, j1);
        float v0  = fy*v01 + (1 - fy)*v00;
        float v1  = fy*v11 + (1 - fy)*v10;
        float v   = fx*v0  + (1 - fx)*v1;
        return v;

    }
};

/// Sample air parcel in the skirt / cap
struct parcel
{
    Eigen::Vector3f position; ///< 3d spatial position
    float vapour;             ///< water vapour value in kg / m^3
    float condensation;       ///< condensed water value in kg / m^3
    float timeSinceExpl;      ///< Time since pressure wave passed by this parcel
    float temp;               ///< Temperature in K, cached for performance
    float pressure;           ///< Pressure used for pressurewave adjustment
};

/// a single thin cloud layer subject to plume uplift
class skirt
{
private:
    // MapFloat * height;           ///< grid of elevation offset values in metres
    // MapFloat * vapour;           ///< grid of water vapour values in kg / m^3
    // MapFloat * condensation;     ///< grid of condensed water values in kg / m^3
    // Eigen::Vector3f start;       ///< starting bottom corner position in metres
    // Eigen::Vector3f extent;      ///< vector extent in metres
    // float step;                  ///< interval between grid vertices in metres
    std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>> atmoModel; ///< access to airlayers
    
    Eigen::Vector3f center;
    float radius;
    float step;
    float thickness;
    float targetheight;
    float condensationThres;
    bool targetreached;
    
    /// Obtain grid size @a dx and @a dy
    // void getGridDim(int & dx, int & dy) const;
    
    /// convert from grid coordinates (gx, gy) to world coordinates (wx, wy, wz)
    // void convertGridToWorld(int gx, int gy, float &wx, float &wy, float &wz);
    
    /// write skirt data to file
    void write(const std::string &filename);
    
    /// initialise skirt with correct water vapour and condensed water from the atmosphere
    // void transferWaterFromAtmosphere(int gx, int gy);

public:
    std::vector<parcel> samples; ///< disk of parcel samples
    
    skirt()
    {
        center = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        radius = 0.0f;
    }

    skirt(Eigen::Vector3f i_center, float i_rad, float i_step, float i_targethght, float i_condenseThres, float i_thickness, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& i_model);
    
    ~skirt()
    {
        samples.clear();
    }
    
    /// initialise layer data
    void init(Eigen::Vector3f i_center, float i_rad, float i_step, float i_targethght, float i_condenseThres, float i_thickness, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& i_model);
    
    /// lift the skirt according to the plume's motion
    /// different areas of influence horizontally (hzone) and vertically (vzone)
    virtual void uplift(scene_model &plume, float dt, float hzone, float vzone);
    
    /// test for condensation against the atmosphere
    void phaseTransition();
    
    /// update temp from atmo
    void updateTemp();

    /// generate file name and call write to save file to disk
    void exportSkirt(int framenum, int layernum);
};

class wilsonCloud : public skirt{
private:
    static constexpr float pressureWaveLife = 1.5;
    float expStart;
    float expLife;

    bool reset, resetOnce;

    std::shared_ptr<PWM::Model::pressureWave<double>> explosion;
public:
    wilsonCloud(Eigen::Vector3f i_center, float i_radius, float i_step, float startTime, float lifeTime, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& i_model);
    ~wilsonCloud();

    void updateExp(float dT);
    //Function to add an explosion to the timeline
    void addExplosion(std::shared_ptr<PWM::Model::pressureWave<double>>& exp);
    void uplift(scene_model &plume, float dt, float hzone, float vzone) override;
    void init(Eigen::Vector3f i_center, float i_rad, float i_step, float i_targethght, float i_condenseThres, float i_thickness, std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& i_model);
};

#endif //CF_SKIRT_HPP
