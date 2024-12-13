#include "atmofuncs.h"
#include "smoke.hpp"

/**
 * Define CONSTVEL to
 */
// #define CONSTVEL
#ifdef CONSTVEL
#define constVelX 10
#define constVelY 0
#endif

// Counter used to save image on hard drive
int counter_image = 0;

// std::chrono::system_clock::time_point a = std::chrono::system_clock::now();
// std::chrono::system_clock::time_point b = std::chrono::system_clock::now();



//------------------------------------------------------------
//------------------------ TERRAIN ---------------------------
//------------------------------------------------------------

void scene_model::fill_height_field(std::vector<Eigen::Vector3f>& position, std::vector<Eigen::Vector3f>& normal,
                                    Eigen::Matrix3f rotation, float scaling, Eigen::Vector3f translation)
{
    std::cerr << "FILL HEIGHT FIELD CALLED" << std::endl;
    // prepare field with parameters
    terrain_struct.min_xyz = -15000.0f;
    terrain_struct.max_xyz = 15000.0f;
    float interval_size = terrain_struct.max_xyz - terrain_struct.min_xyz;
    terrain_struct.cell_size = 200.f;
    terrain_struct.field_size = (size_t)(terrain_struct.max_xyz/terrain_struct.cell_size);
    terrain_struct.height_field.resize(terrain_struct.field_size);
    terrain_struct.normal_field.resize(terrain_struct.field_size);
    for (unsigned int i = 0; i < terrain_struct.field_size; i++)
    {
        terrain_struct.height_field[i].resize(terrain_struct.field_size);
        terrain_struct.normal_field[i].resize(terrain_struct.field_size);
    }

    //transform like for mesh_drawable
    for (unsigned int i = 0; i<position.size(); i++)
    {
        position[i] = scaling * (rotation * position[i] + translation);
        normal[i] = rotation * normal[i];
    }
    terrain_struct.positions = position;
    terrain_struct.normals = normal;


    //fill height field for collisions
    for (unsigned int i = 0; i<position.size(); i++)
    {
        if (position[i].x() < terrain_struct.max_xyz && position[i].x() > terrain_struct.min_xyz
                && position[i].y() < terrain_struct.max_xyz && position[i].y() > terrain_struct.min_xyz)
        {
            int idx_x = (int)(terrain_struct.field_size*(position[i].x() - terrain_struct.min_xyz)/interval_size);
            int idx_y = (int)(terrain_struct.field_size*(position[i].y() - terrain_struct.min_xyz)/interval_size);
            if (idx_x == terrain_struct.field_size) idx_x = terrain_struct.field_size-1;
            if (idx_y == terrain_struct.field_size) idx_y = terrain_struct.field_size-1;
            terrain_struct.height_field[idx_x][idx_y] = position[i].z();
            terrain_struct.normal_field[idx_x][idx_y] = normal[i];
        }
    }

}

void scene_model::addAtmoModel(std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& model){
    this->atmoModel = model;
}

void scene_model::setPlumeLoc(Eigen::Vector3f loc){
    this->plume_loc = loc;
    this->z_0 = loc.z();
}

//Example use of world model to sample velocity:
Eigen::Vector2f scene_model::sampleVel(Eigen::Vector3f& loc){
    #ifdef CONSTVEL
        return Eigen::Vector2f(constVelX, constVelY);
    #endif
    auto lVel = atmoModel->getLocalVel(loc.x(), loc.y(), loc.z());
    return Eigen::Vector2f(lVel.first, lVel.second);
}

/// Skirt cloud helper functions


float scene_model::calcFallOff(Eigen::Vector3f& loc, Eigen::Vector3f& c, float r, float hzone, float vzone)
{
    float t, dist, wght = 0.0f;
    Eigen::Vector3f cloc = loc - c;
    
    dist = cloc.squaredNorm();
    float rsq = r * r;
    if(dist <= rsq)
    {
        wght = 1.0f;
    }
    else
    {
        // ellipsoidal area of influence
        Eigen::Vector3f dvec = Eigen::Vector3f(cloc.x() / hzone, cloc.y() / hzone, cloc.z() / vzone);
        dist = dvec.squaredNorm();
        if(dist < 1.0f)
        {
            // smoothstep for decay
            // float t = 1.0 - (dist - r) / r;
            // float t = 1.0 - (dist - rsq)/ zonesq;
            dist = sqrt(dist);
            t = 1.0 - dist;
            wght = t * t * (3.0f - 2.0f * t);
        }
    }
   
    /*
    float zonesq = zone * zone;
    if(dist <= rsq + zonesq)
    {
        if(dist <= rsq)
        {
            wght = 1.0f;
        }
        else
        {
            // smoothstep for decay
            // float t = 1.0 - (dist - r) / r;
            float t = 1.0 - (dist - rsq)/ zonesq;
            wght = t * t * (3.0f - 2.0f * t);
        }
    }*/
    
    return wght;
}


void scene_model::boxSpheres()
{
    float minx = 60000.0f, miny = 60000.0f, maxx = -60000.0f, maxy = -60000.0f, minz = 60000.0f, maxz = -60000.0f;
    
    for (int i = 0; i < smoke_layers.size(); i++)
    {
        if(smoke_layers[i].center.x() < minx)
            minx = smoke_layers[i].center.x();
        if(smoke_layers[i].center.x() > maxx)
            maxx = smoke_layers[i].center.x();
        if(smoke_layers[i].center.y() < miny)
            miny = smoke_layers[i].center.y();
        if(smoke_layers[i].center.y() > maxy)
            maxy = smoke_layers[i].center.y();
        if(smoke_layers[i].center.z() < minz)
            minz = smoke_layers[i].center.z();
        if(smoke_layers[i].center.z() > maxz)
            maxz = smoke_layers[i].center.z();
    }
    
    for (int i = 0; i < free_spheres.size(); i++)
    {
        if(free_spheres[i].center.x() < minx)
            minx = free_spheres[i].center.x();
        if(free_spheres[i].center.x() > maxx)
            maxx = free_spheres[i].center.x();
        if(free_spheres[i].center.y() < miny)
            miny = free_spheres[i].center.y();
        if(free_spheres[i].center.y() > maxy)
            maxy = free_spheres[i].center.y();
        if(free_spheres[i].center.z() < minz)
            minz = free_spheres[i].center.z();
        if(free_spheres[i].center.z() > maxz)
            maxz = free_spheres[i].center.z();
       
    }
    
    for (int i = 0; i < s2_spheres.size(); i++)
    {
        if(s2_spheres[i].center.x() < minx)
            minx = s2_spheres[i].center.x();
        if(s2_spheres[i].center.x() > maxx)
            maxx = s2_spheres[i].center.x();
        if(s2_spheres[i].center.y() < miny)
            miny = s2_spheres[i].center.y();
        if(s2_spheres[i].center.y() > maxy)
            maxy = s2_spheres[i].center.y();
        if(s2_spheres[i].center.z() < minz)
            minz = s2_spheres[i].center.z();
        if(s2_spheres[i].center.z() > maxz)
            maxz = s2_spheres[i].center.z();
       
    }
    
    for (int i = 0; i < s3_spheres.size(); i++)
    {
        if(smoke_layers[i].center.x() < minx)
            minx = smoke_layers[i].center.x();
        if(smoke_layers[i].center.x() > maxx)
            maxx = smoke_layers[i].center.x();
        if(smoke_layers[i].center.y() < miny)
            miny = smoke_layers[i].center.y();
        if(smoke_layers[i].center.y() > maxy)
            maxy = smoke_layers[i].center.y();
        if(smoke_layers[i].center.z() < minz)
            minz = smoke_layers[i].center.z();
        if(smoke_layers[i].center.z() > maxz)
            maxz = smoke_layers[i].center.z();
       
    }
    
    std::cerr << "PLUME BOUNDS" << std::endl;
    std::cerr << "X [" << minx << ", " << maxx << "]" << std::endl;
    std::cerr << "Y [" << miny << ", " << maxy << "]" << std::endl;
    std::cerr << "Z [" << minz << ", " << maxz << "]" << std::endl;
}

void scene_model::boundSpheres(float zone)
{
    for (int i = 0; i < smoke_layers.size(); i++)
    {
        smoke_layers[i].maxz = smoke_layers[i].center.z() + smoke_layers[i].r + zone;
        smoke_layers[i].minz = smoke_layers[i].center.z() - smoke_layers[i].v.z() - smoke_layers[i].r - zone;
    }
    
    for (int i = 0; i < free_spheres.size(); i++)
    {
        free_spheres[i].maxz = free_spheres[i].center.z() + free_spheres[i].r + zone;
        free_spheres[i].minz = free_spheres[i].center.z() - free_spheres[i].speed.z() - free_spheres[i].r - zone;
    }
    
    for (int i = 0; i < s2_spheres.size(); i++)
    {
        s2_spheres[i].maxz = s2_spheres[i].center.z() + s2_spheres[i].r + zone;
        s2_spheres[i].minz = s2_spheres[i].center.z() - s2_spheres[i].speed.z() - s2_spheres[i].r - zone;
    }
    
    for (int i = 0; i < s3_spheres.size(); i++)
    {
        s3_spheres[i].maxz = s3_spheres[i].center.z() + s3_spheres[i].r + zone;
        s3_spheres[i].minz = s3_spheres[i].center.z() - s3_spheres[i].speed.z() - s3_spheres[i].r - zone;
    }
}

int scene_model::calcSubsteps()
{
    int steps, maxsteps = 1;
    
    for (int i = 0; i < smoke_layers.size(); i++)
    {
        steps = (int) (smoke_layers[i].v.norm() / (2.0f * smoke_layers[i].r));
        if(steps > maxsteps)
        {
            if(steps > 500)
                std::cerr << "Too many substeps"<< std::endl;
            maxsteps = steps;
        }
    }
    for (int i = 0; i < free_spheres.size(); i++)
    {
        steps = (int) (free_spheres[i].speed.norm() / (2.0f * free_spheres[i].r));
        if(steps > maxsteps)
        {
            if(steps > 500)
                std::cerr << "Too many substeps"<< std::endl;
            maxsteps = steps;
        }
    }
    /*
    for (int i = 0; i < s2_spheres.size(); i++)
    {
        steps = (int) (s2_spheres[i].speed.norm() / (2.0f * s2_spheres[i].r));
        if(steps > maxsteps)
        {
            if(steps > 500)
                std::cerr << "Too many substeps"<< std::endl;
            maxsteps = steps;
        }
    }*/
    return maxsteps;
}

void scene_model::trajectoryInterp(Eigen::Vector3f& loc, Eigen::Vector3f& vel, Eigen::Vector3f& subloc, Eigen::Vector3f& subvel, float t)
{
    Eigen::Vector3f start;
    
    subvel = t * vel;
    start = loc - vel;
    subloc = start + subvel;
}

float scene_model::calcSkirtUplift(Eigen::Vector3f& loc, float deltaT, float hzone, float vzone, int substeps)
{
    float lift, totlift = 0.0f;
    Eigen::Vector3f subcenter, subvel;
    
    // find maximum lift among all slices and attached spheres, subdivided into substeps
    
    float tstep = 1.0f / (float) (substeps+1);
    for(int s = 0; s <= substeps; s++)
    {
        float t, maxlift = 0.0f;
        
        // interpolation position
        if(s == substeps)
            t = 1.0f;
        else
            t = tstep * (float) (s+1);
        
        for (int i = 0; i < smoke_layers.size(); i++)
        {
            if(loc.z() >= smoke_layers[i].minz && loc.z() <= smoke_layers[i].maxz)
            {
                trajectoryInterp(smoke_layers[i].center, smoke_layers[i].v, subcenter, subvel, t);
                lift = subvel.z() * calcFallOff(loc, subcenter, smoke_layers[i].r, hzone, vzone);
                if(lift > maxlift)
                    maxlift = lift;
            }
        }
        
        for (int i = 0; i < free_spheres.size(); i++)
        {
            if(loc.z() >= free_spheres[i].minz && loc.z() <= free_spheres[i].maxz)
            {
                trajectoryInterp(free_spheres[i].center, free_spheres[i].speed, subcenter, subvel, t);
                lift = subvel.z() * deltaT * calcFallOff(loc, subcenter, free_spheres[i].r, hzone, vzone);
                if(lift > maxlift)
                    maxlift = lift;
            }
        }
        
        
        for (int i = 0; i < s2_spheres.size(); i++)
        {
            if(loc.z() >= s2_spheres[i].minz && loc.z() <= s2_spheres[i].maxz)
            {
                trajectoryInterp(s2_spheres[i].center, s2_spheres[i].speed, subcenter, subvel, t);
                lift = subvel.z() * deltaT * calcFallOff(loc, subcenter, s2_spheres[i].r, hzone, vzone);
                if(lift > maxlift)
                    maxlift = lift;
            }
        }
/*
        for (int i = 0; i < s3_spheres.size(); i++)
        {
            if(loc.z() >= s2_spheres[i].minz && loc.z() <= s2_spheres[i].maxz)
            {
                trajectoryInterp(s3_spheres[i].center, s3_spheres[i].speed, subcenter, subvel, t);
                lift = subvel.z() * deltaT * calcFallOff(loc, subcenter, s3_spheres[i].r, hzone, vzone);
                if(lift > maxlift)
                    maxlift = lift;
            }
        }*/
        totlift += maxlift;
    }
    
    return totlift;
}

float scene_model::truncNorm(float val, float max)
{
    return fminf(max, val) / max;
}
 
// calculate the change in delta for a single time-step
float scene_model::chargeDelta(Eigen::Vector3f& loc, Eigen::Vector3f& vel, float density)
{
    float delta;
    
    // collect attributes that affect charge
 
    // atmospheric temperature
    float atmp = atmoModel->getLocalTemp(loc.x(), loc.y(), loc.z());
    
    if(atmp >= 253.15f) // behaviour depends on -20C isotherm
    {
        // combination of ash and sphere velocity, suppressed by water content
 
        float ash = truncNorm(density, ashChargeHotMax);
        float speed = truncNorm(vel.norm(), velChargeMax);
        float moisture = atmoModel->getLocalMoist(loc.x(), loc.y(), loc.z());
        float moistwght = exp(-10.0f * truncNorm(moisture, moistureChargeMax));
        
        delta = hotChargeCoeff * ash * speed * moistwght;
        chrghotcnt++;
        chrghotsum += delta;
      
        return delta;
        
    }
    else // ice starts to play a significant role
    {
        // combination of ash, water(ice), convection
 
        float ice = atmoModel->getLocalClouds(loc.x(), loc.y(), loc.z());
        ice = truncNorm(ice, iceChargeMax);
        float ash = truncNorm(density, ashChargeColdMax);

        float advect = fabs(atmoModel->getLocalVelX(loc.x(), loc.y(), loc.z()))+fabs(atmoModel->getLocalVelY(loc.x(), loc.y(), loc.z()));
        #ifdef CONSTVEL
            advect = std::abs(constVelX) + std::abs(constVelY);
        #endif
        float convect = fabs(atmoModel->getLocalConvection(loc.x(), loc.y(), loc.z()));
        convect = truncNorm(convect+advect, convectChargeMax);
        
        delta = coldChargeCoeff * convect * 0.5f * (ash + ice);
        chrgcoldcnt++;
        chrgcoldsum += delta;
        
        return delta;
    }
}

//------------------------------------------------------------
//------------------------ PROCESS ---------------------------
//------------------------------------------------------------ */

void scene_model::frame_draw(float deltaT)
{
    // Maintain designated frequency of 5 Hz (200 ms per frame)
//    a = std::chrono::system_clock::now();
//    std::chrono::duration<double, std::milli> work_time = a - b;

//    if (work_time.count() < 16.0)
//    {
//        std::chrono::duration<double, std::milli> delta_ms(16.0 - work_time.count());
//        auto delta_ms_duration = std::chrono::duration_cast<std::chrono::milliseconds>(delta_ms);
//        std::this_thread::sleep_for(std::chrono::milliseconds(delta_ms_duration.count()));
//    }

//    b = std::chrono::system_clock::now();
//    std::chrono::duration<double, std::milli> sleep_time = b - a;

    // const float dt = timer.update();
    //set_gui();

    // Force constant time step
    // dt = dt<=1e-6f? 0.0f : timer.scale*0.002f; //0.0003f
    auto a = std::chrono::system_clock::now();
    new_layer_delay += deltaT;
    dt = deltaT;

    // add smoke layer each x seconds
    if (smoke_layers.size() == 0 || (new_layer_delay >= r_0/(2*U_0) && smoke_layers.size() < 1000000000000000))
    {
        add_smoke_layer(U_0, rho_0, c_0, r_0, plume_loc, false);
        add_free_spheres_for_one_layer(smoke_layers.size()-1);

        new_layer_delay = 0;
        total_layers_ejected++;
        //std::cout << smoke_layers.size() << std::endl;
        //if (debug_mode) std::cout << "LAYER ADDED OK" << std::endl;
    }

    chrghotsum = 0.0f;
    chrghotcnt = 0;
    chrgcoldsum = 0.0f;
    chrgcoldcnt = 0;
    
    // update of layer and spheres
//#pragma omp parallel for
    for (unsigned int id=0; id<smoke_layers.size(); id++)
    {
        smoke_layer_update(id);
    }
//#pragma omp barrier
    
    update_free_spheres();
    if (frame_count %100 == 0) falling_spheres_update(100);
    update_stagnation_spheres();
    
    update_sphere_charges();
    
    /*
    if(frame_count % 200 == 0)
    {
        if(chrghotcnt > 0)
            std::cerr << "hot avg = " << chrghotsum / (float) chrghotcnt << ", count = " << chrghotcnt << std::endl;
        if(chrgcoldcnt > 0)
            std::cerr << "cold avg = " << chrgcoldsum / (float) chrgcoldcnt << ", count = " << chrgcoldcnt << std::endl;
    }*/

    // update subspheres
    //if (frame_count %50 == 0) update_subspheres_params();

    // export (comment or uncomment)
 //   if (export_data && frame_count % 50 == 0)
 //   {
 //       export_spheres_combined(50);
//        export_spheres(50);
//        export_spheres_binary(50);
 //   }
    /*
    // store data for replay
    if (!export_data && frame_count %50 == 0)
    {
        smoke_layers_frames.push_back(smoke_layers);
        free_spheres_frames.push_back(free_spheres);
        stagnate_spheres_frames.push_back(stagnate_spheres);
        falling_spheres_frames.push_back(falling_spheres);
        falling_spheres_buffers_frames.push_back(falling_spheres_buffers);
    }*/

    frame_count++;
    t_sum += dt;

    auto b = std::chrono::system_clock::now();
    run_time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count()/1000000000.0;
    /*
    if (replay)
    {
        display_replay(shaders, scene, gui, frame_replay);
        std::cout << frame_replay << std::endl;
        frame_replay++;
        frame_replay = frame_replay % smoke_layers_frames.size();
    }
    else display(shaders, scene, gui);*/
    
    // report terrain statistics
    /*
    std::cerr << "TERRAIN STATS" << std::endl;
    std::cerr << "minxyz = " << terrain_struct.min_xyz << " maxxyz = " << terrain_struct.max_xyz << std::endl;
     */

}


//------------------------------------------------------------
//------------------------- ALGO -----------------------------
//------------------------------------------------------------ */

void scene_model::add_smoke_layer(float v, float d, float c, float r, Eigen::Vector3f position, bool secondary_plume)
{
    smoke_layer layer = smoke_layer({0,0,v}, d, c, r, position, secondary_plume);
    layer.moisture = 10.0f; // JG - put water vapour at the vent point, 0.5
    smoke_layers.push_back(layer);
}

float scene_model::compute_atm_temperature(float height) //(K) cf https://protect-za.mimecast.com/s/H4ByCMjKp1Sq2JDXlC3kmNn
{
    return 288.15 - 6.5*height/1000.0;
}

float scene_model::compute_atm_density(float height) //(kg.m-3) cf https://protect-za.mimecast.com/s/UUB1CNxKq0i0Z6J7RcpdKPs
{
    return 352.995 * pow(1 - 0.0000225577*height, 5.25516) / (288.15 - 0.0065*height);
}

float scene_model::compute_atm_pressure(float height) //(Pa) cf https://protect-za.mimecast.com/s/Eaw7CO7Xr8spNOlnGSNu5rB
{
    return 101325. * pow(1. - 0.0065 * height / 288.15, 5.255);
}

Eigen::Vector3f scene_model::compute_wind_vector(float height)
{
    // find altitude interval
    unsigned int low_altitude_idx = 0;
    for (unsigned int i = 0; i<wind_altitudes.size(); i++)
    {
        if (wind_altitudes[i] < height)
        {
            low_altitude_idx = i;
        }
    }

    // compute wind vec by interpolating
    if (low_altitude_idx == wind_altitudes.size()-1)
    {
        return winds[low_altitude_idx].wind_vector;
    }
    else
    {
        float low_height = (float)wind_altitudes[low_altitude_idx];
        float high_height = (float)wind_altitudes[low_altitude_idx+1];
        float lambda = (height-low_height)/(high_height-low_height);
        Eigen::Vector3f interpo_wind = winds[low_altitude_idx].wind_vector + lambda * (winds[low_altitude_idx+1].wind_vector - winds[low_altitude_idx].wind_vector);
        return interpo_wind;
    }
}

float scene_model::compute_temperature_at_h(float height)
{
    // find altitude interval
    unsigned int low_altitude_idx = 0;
    for (unsigned int i = 0; i<wind_altitudes.size(); i++)
    {
        if (wind_altitudes[i] < height)
        {
            low_altitude_idx = i;
        }
    }

    // compute wind vec by interpolating
    if (low_altitude_idx == wind_altitudes.size()-1)
    {
        return temperatures[low_altitude_idx];
    }
    else
    {
        float low_height = (float)wind_altitudes[low_altitude_idx];
        float high_height = (float)wind_altitudes[low_altitude_idx+1];
        float lambda = (height-low_height)/(high_height-low_height);
        float interpo_temp = temperatures[low_altitude_idx] + lambda * (temperatures[low_altitude_idx+1] - temperatures[low_altitude_idx]);
        return interpo_temp;
    }
}

float scene_model::compute_moisture_at_h(float height)
{
    // find altitude interval
    unsigned int low_altitude_idx = 0;
    for (unsigned int i = 0; i<wind_altitudes.size(); i++)
    {
        if (wind_altitudes[i] < height)
        {
            low_altitude_idx = i;
        }
    }

    // compute wind vec by interpolating
    if (low_altitude_idx == wind_altitudes.size()-1)
    {
        return moistures[low_altitude_idx];
    }
    else
    {
        float low_height = (float)wind_altitudes[low_altitude_idx];
        float high_height = (float)wind_altitudes[low_altitude_idx+1];
        float lambda = (height-low_height)/(high_height-low_height);
        float interpo_temp = moistures[low_altitude_idx] + lambda * (moistures[low_altitude_idx+1] - moistures[low_altitude_idx]);
        return interpo_temp;
    }
}

float scene_model::compute_density_at_h(float height) // rho = M/R * P/T
{
    // find altitude interval
    unsigned int low_altitude_idx = 0;
    for (unsigned int i = 0; i<wind_altitudes.size(); i++)
    {
        if (wind_altitudes[i] < height)
        {
            low_altitude_idx = i;
        }
    }

    // compute wind vec by interpolating
    if (low_altitude_idx == wind_altitudes.size()-1)
    {
        return PWM::Utils::pressureAdjustedDensity(pressures[low_altitude_idx], temperatures[low_altitude_idx], planetos);
    }
    else
    {
        float low_height = (float)wind_altitudes[low_altitude_idx];
        float high_height = (float)wind_altitudes[low_altitude_idx+1];
        float lambda = (height-low_height)/(high_height-low_height);
        float interpo_pres = pressures[low_altitude_idx] + lambda * (pressures[low_altitude_idx+1] - pressures[low_altitude_idx]);
        float interpo_temp = temperatures[low_altitude_idx] + lambda * (temperatures[low_altitude_idx+1] - temperatures[low_altitude_idx]);
        return PWM::Utils::pressureAdjustedDensity(interpo_pres, interpo_temp, planetos);
    }
}

float scene_model::compute_pressure_at_h(float height)
{
    // find altitude interval
    unsigned int low_altitude_idx = 0;
    for (unsigned int i = 0; i<wind_altitudes.size(); i++)
    {
        if (wind_altitudes[i] < height)
        {
            low_altitude_idx = i;
        }
    }

    // compute wind vec by interpolating
    if (low_altitude_idx == wind_altitudes.size()-1)
    {
        return pressures[low_altitude_idx];
    }
    else
    {
        float low_height = (float)wind_altitudes[low_altitude_idx];
        float high_height = (float)wind_altitudes[low_altitude_idx+1];
        float lambda = (height-low_height)/(high_height-low_height);
        float interpo_pres = pressures[low_altitude_idx] + lambda * (pressures[low_altitude_idx+1] - pressures[low_altitude_idx]);
        return interpo_pres;
    }
}

void scene_model::edit_smoke_layer_properties(unsigned int i, float& d_mass)
{
    // get wind at altitude
    Eigen::Vector3f wind = compute_wind_vector(smoke_layers[i].center.z());

    // preliminary computation
    float thk = dt * smoke_layers[i].v.z();
    //thk = smoke_layers[i].v.z * 0.002;

    float total_smoke_thk = smoke_layers[i].thickness;
    float total_smoke_volume = total_smoke_thk * 3.14 * smoke_layers[i].r * smoke_layers[i].r;
    float total_smoke_mass = smoke_layers[i].rho * total_smoke_volume;

    // wind velocity around for air incorporation
    float k_s = 0.09, k_w = 0.9;
    float U_e = k_s * abs((smoke_layers[i].v).norm() - wind.norm() * cos(smoke_layers[i].theta)) + k_w * abs(wind.norm() * sin(smoke_layers[i].theta));
    float r_atm = air_incorporation_coeff * U_e;
    //r_atm = U_e * dt;

    // air quantity to put in the plume
    float atm_density = compute_atm_density(smoke_layers[i].center.z());
    //float atm_density = compute_density_at_h(smoke_layers[i].center.z());
    float atm_volume = thk * 3.14 * (2.0 * smoke_layers[i].r + r_atm) * r_atm; //air around the smoke layer
    float atm_mass = atm_density * atm_volume;

    // compute new temperature (we take all Cp equal)
    //float atm_temperature = compute_atm_temperature(smoke_layers[i].center.z());
    float atm_temperature = compute_temperature_at_h(smoke_layers[i].center.z());
    float new_temp = (total_smoke_mass * smoke_layers[i].temperature + atm_mass * atm_temperature) / (total_smoke_mass + atm_mass);
    smoke_layers[i].temperature = new_temp;

    // new volume after heating
    float atm_new_volume = smoke_layers[i].temperature * atm_volume / atm_temperature; // after heating by hot smoke

    // moisture
    float water_quantity = compute_moisture_at_h(smoke_layers[i].center.z())*atm_volume;

    // new params
    d_mass = atm_mass;
    float mass_new = total_smoke_mass + atm_mass;
    float volume_new = total_smoke_volume + atm_new_volume;
    float rho_new = mass_new/volume_new;
    float moist_new = water_quantity/volume_new;
    float r_new = cbrt(volume_new/3.14);
    smoke_layers[i].thickness = r_new;


    // new speed due to conservation of energy (old)
    //float energy = 0.5 * total_smoke_mass * smoke_layers[i].v.z * smoke_layers[i].v.z;
    //float new_speed = sqrt(2.0 * energy / mass_new);
    //if (rho_new < atm_density) smoke_layers[i].v = {0,0,new_speed};

    if (smoke_layers[i].rho > atm_density && rho_new < atm_density) smoke_layers[i].plume = true;
    smoke_layers[i].r = r_new;
    smoke_layers[i].rho = rho_new;
    smoke_layers[i].moisture += moist_new;
}

void scene_model::apply_forces_to_smoke_layer(unsigned int i, float d_mass)
{
    // get wind at altitude
    Eigen::Vector3f wind = compute_wind_vector(smoke_layers[i].center.z());

    // precomputation
    float atm_density = compute_atm_density(smoke_layers[i].center.z());
    //float atm_density = compute_density_at_h(smoke_layers[i].center.z());
    float volume = smoke_layers[i].thickness*3.14*smoke_layers[i].r*smoke_layers[i].r;
    float surface = 2*3.14*smoke_layers[i].r*smoke_layers[i].r;
    float surface_eff = 2*smoke_layers[i].r*smoke_layers[i].r;
    float smoke_mass = smoke_layers[i].rho * volume;
    Eigen::Vector3f v_diff = wind - Eigen::Vector3f(smoke_layers[i].v.x(), smoke_layers[i].v.y(), 0);

    Eigen::Vector3f weight = {0,0, -smoke_mass * g}; // m*g
    Eigen::Vector3f archimede = {0,0, atm_density * volume * g}; // rho*V*g
    Eigen::Vector3f friction = - 0.5 * atm_density * 0.04 * surface * (smoke_layers[i].v).norm() * smoke_layers[i].v; // axial friction
    //friction = - volume * 0.00005 * norm(smoke_layers[i].v) * smoke_layers[i].v; // old way to compute friction
    Eigen::Vector3f wind_force = 600. * (v_diff).norm() * v_diff * surface_eff; // horizontal

    Eigen::Vector3f forces = weight + archimede + friction + wind_force;

    Eigen::Vector3f old_v = smoke_layers[i].v;
    float old_m = smoke_mass - d_mass;

    // equation of dynamics without mass conservation
    Eigen::Vector3f mv = old_m * old_v + forces*dt;
    Eigen::Vector3f v = mv/smoke_mass;
    Eigen::Vector3f p = smoke_layers[i].center + v*dt;

    // old with mass conservation
    //vec3 a = forces / (smoke_mass);
    //vec3 v = smoke_layers[i].v + a*dt;
    //vec3 p = smoke_layers[i].center + v*dt;


    // check if falls
    if (old_v.z() >0 && v.z() < 0 && smoke_layers[i].plume == false)
    {
        smoke_layers[i].rising = false;
        smoke_layers[i].begin_falling = true;
        std::cout << "layer " << i << " falling frame " << frame_count << std::endl;
    }

    // check if stagnates
    if (smoke_layers[i].plume && !smoke_layers[i].stagnates && smoke_layers[i].center.z() > 5000 && (atm_density - smoke_layers[i].rho < 0.01))
    {
        smoke_layers[i].stagnates = true;
        stagnation_speed = smoke_layers[i].v.z();
    }

    // check if stagnates long
    if (smoke_layers[i].stagnates && v.z() < 0 && !smoke_layers[i].stagnates_long)
    {
        smoke_layers[i].stagnates_long = true;
    }

    if (smoke_layers[i].stagnates && p.z() < smoke_layers[i].center.z())
    {
        p.z() = smoke_layers[i].center.z();
    }

    // update
    smoke_layers[i].a = v/dt;
    smoke_layers[i].v = v;
    smoke_layers[i].center = p;
    if (p.z() < -1000) smoke_layers[i].center.z() = -1000; // prevent from going oob after falling (layers are not deleted but not used anymore)

    // compute new theta
    // WARNING: wind direction was constant in my tests, if the direction changes it may not work, it has to be tested
    //float dx = v.x*dt, dy = v.y*dt;
    float dz = v.z() * dt;
    smoke_layers[i].speed_along_axis = v.norm();
    smoke_layers[i].plume_axis = v.normalized();
    if (wind.norm() != 0 && dt > 0)
    {
        smoke_layers[i].theta_axis = (wind.cross(Eigen::Vector3f(0,0,1))).normalized();
        float theta_totest = asin(dz/(v*dt).norm());
        if (dz/(v*dt).norm() < 1.000001 && dz / (v * dt).norm() > 0.999999) theta_totest = 3.14 / 2.; // security to prevent nan values due to asin
        smoke_layers[i].theta = theta_totest;

        if (theta_totest < 0)
        {
            smoke_layers[i].theta = 0;
        }
    }


    // print values in txt files
    nb_of_iterations++;
    if (i == 0)
    {
        if (altitude) altitude << smoke_layers[i].center.z() << std::endl;
        if (plumex) plumex << smoke_layers[i].center.x() << std::endl;
        if (speed) speed << smoke_layers[i].v.z() << std::endl;
        if (ray) ray << smoke_layers[i].r << std::endl;
        if (smoke_rho) smoke_rho << smoke_layers[i].rho << std::endl;
        if (atm_rho) atm_rho << atm_density << std::endl;
        if (temp) temp << smoke_layers[i].temperature << std::endl;
        if (moist) moist << smoke_layers[i].moisture << std::endl;
        if (time) time << t_sum << std::endl;

        if (theoric_temp) theoric_temp << compute_atm_temperature(smoke_layers[i].center.z()) << std::endl;
        //if (actual_temp) actual_temp << compute_temperature_at_h(smoke_layers[i].center.z()) << std::endl;
        if (actual_temp) actual_temp << compute_temperature_at_h(2500) << std::endl;
        /*if (theoric_rho) theoric_rho << compute_atm_density(smoke_layers[i].center.z()) << std::endl;
        if (actual_rho) actual_rho << compute_density_at_h(smoke_layers[i].center.z()) << std::endl;
        if (theoric_pres) theoric_pres << compute_atm_pressure(smoke_layers[i].center.z()) << std::endl;
        if (actual_pres) actual_pres << compute_pressure_at_h(smoke_layers[i].center.z()) << std::endl;*/
    }
    
    // JG v * dt == change in p?

}

void scene_model::sedimentation(unsigned int i, float& d_mass)
{
    // constant sedimentation
    float layer_volume = 3.14 * smoke_layers[i].r*smoke_layers[i].r * smoke_layers[i].thickness;
    float diff_density = 0.00000005 * dt;
    if (smoke_layers[i].stagnates_long) diff_density = 0.00005 * dt;

    if (smoke_layers[i].rho > diff_density)
    {
        smoke_layers[i].rho -= diff_density;
        d_mass -= diff_density * layer_volume;
    }

}

void scene_model::smoke_layer_update(unsigned int i)
{
    float d_mass = 0; // to track mass change for equation of dynamics
    if (smoke_layers[i].plume == true && smoke_layers[i].center.z() > 0.) sedimentation(i, d_mass); // sedimentation in altitude
    if (smoke_layers[i].rising && !smoke_layers[i].stagnates_long) edit_smoke_layer_properties(i, d_mass); // convection if v_z > 0 (convection causes air entrainment)
    apply_forces_to_smoke_layer(i, d_mass);
    
    float delta = chargeDelta(smoke_layers[i].center, smoke_layers[i].v, smoke_layers[i].rho);
    smoke_layers[i].charge += delta; // increment charge
}

void scene_model::update_sphere_charges()
{
    int freeSpheresSize = free_spheres.size() - 1;
#pragma omp parallel for
    for (int i = freeSpheresSize; i>=0; i--)
    {
        free_spheres[i].charge += chargeDelta(free_spheres[i].center, free_spheres[i].speed, free_spheres[i].rho);
    }
#pragma omp barrier

#pragma omp parallel for
    // borrow density from parent sphere
    for (int i = 0; i<s2_spheres.size(); i++)
    {
        float parentrho = free_spheres[s2_spheres[i].parent_id].rho;
        s2_spheres[i].charge += chargeDelta(s2_spheres[i].center, s2_spheres[i].speed, parentrho);
    }
#pragma omp barrier

    int fallSpheresSize = falling_spheres.size() - 1;
#pragma omp parallel for
    for (int i = fallSpheresSize; i>=0; i--)
    {
        falling_spheres[i].charge += chargeDelta(falling_spheres[i].center, falling_spheres[i].speed, falling_spheres[i].rho);
        
    }
#pragma omp barrier
}

//------------------------------------------------------------
//--------------------- FALLING SPHERES ----------------------
//------------------------------------------------------------ */

float scene_model::field_height_at(float x, float y)
{
    float interval_size = terrain_struct.max_xyz - terrain_struct.min_xyz;
    int idx_x = (int)(terrain_struct.field_size*(x - terrain_struct.min_xyz)/interval_size);
    int idx_y = (int)(terrain_struct.field_size*(y - terrain_struct.min_xyz)/interval_size);
    return terrain_struct.height_field[idx_x][idx_y];
}

Eigen::Vector3f scene_model::field_normal_at(float x, float y)
{
    float interval_size = terrain_struct.max_xyz - terrain_struct.min_xyz;
    int idx_x = (int)(terrain_struct.field_size*(x - terrain_struct.min_xyz)/interval_size);
    int idx_y = (int)(terrain_struct.field_size*(y - terrain_struct.min_xyz)/interval_size);
    return terrain_struct.normal_field[idx_x][idx_y];
}

void scene_model::sphere_ground_collision(free_sphere_params& sphere, int idx, unsigned int frame_nb)
{
    // find ground point and normal
    float terrain_z = field_height_at(sphere.center.x(), sphere.center.y());

    // if sphere under ground mesh
    if (sphere.center.z() < terrain_z)
    {

        Eigen::Vector3f terrain_pt(sphere.center.x(), sphere.center.y(), terrain_z);
        Eigen::Vector3f terrain_normal = field_normal_at(sphere.center.x(), sphere.center.y());
        Eigen::Vector3f pt_diff = terrain_pt - sphere.center;

        // compute new position
        sphere.center += pt_diff.norm() * terrain_normal;

        // compute new speed
        Eigen::Vector3f v_normal = (sphere.speed).dot(terrain_normal)* terrain_normal;
        Eigen::Vector3f v_tan = sphere.speed - v_normal;
        sphere.speed = 1.0f*v_tan - 0.0f*v_normal;

        // sedimentation
        if (sphere.falling_under_atm_rho == false) sphere.rho -= 0.001*frame_nb*dt *(sphere.speed).norm();
        else sphere.rho -= 0.000005*frame_nb*dt * (sphere.speed).norm();
        // find buffer for low density particles not in buffer:
        float atm_density = compute_atm_density(sphere.center.z());
        //float atm_density = compute_density_at_h(sphere.center.z());
        if (idx >=0 && sphere.rho < atm_density)
        {
            sphere.falling_under_atm_rho = true;

            // add particle in a buffer

            // test proximity with existing buffers
            bool is_in_buffer = false;
            for (unsigned int j = 0; j<falling_spheres_buffers.size(); j++)
            {
                if ((falling_spheres_buffers[j][0].center - sphere.center).norm() < 5 * sphere.r)
                {
                    falling_spheres_buffers[j].push_back(sphere);
                    falling_spheres.erase(falling_spheres.begin()+idx);
                    is_in_buffer = true;
                    break;
                }
            }
            // if not cloase to any existing buffer, create a new one
            if (is_in_buffer == false)
            {
                std::vector<free_sphere_params> new_buffer;
                new_buffer.push_back(sphere);
                falling_spheres_buffers.push_back(new_buffer);
                falling_spheres.erase(falling_spheres.begin()+idx);
            }
        }
    }

}

void scene_model::ground_falling_sphere_update(free_sphere_params& sphere, int idx, unsigned int frame_nb)
{
    // find closest layer for radial force
    unsigned int closest_layer_id = 0;
    float min_dist = (sphere.center - smoke_layers[0].center).norm();
    for (unsigned int j = 0; j<smoke_layers.size(); j++)
    {
        float dist = (sphere.center - smoke_layers[j].center).norm();
        if (dist < min_dist && smoke_layers[j].rising && smoke_layers[j].secondary_plume == false)
        {
            min_dist = dist;
            closest_layer_id = j;
        }
    }

    // precomputation
    float V = 4./3. * 3.14 * sphere.r * sphere.r * sphere.r;
    float m = sphere.rho * V;

    float atm_rho = compute_atm_density(sphere.center.z());
    //float atm_rho = compute_density_at_h(sphere.center.z());
    Eigen::Vector3f gravity = Eigen::Vector3f(0,0,-m*g);
    Eigen::Vector3f buoyancy = Eigen::Vector3f(0,0,atm_rho*V*g);
    Eigen::Vector3f friction = - 0.1 * (sphere.speed).normalized();
    friction = Eigen::Vector3f(0,0,0);

    // tests for additional force to prevent from being to close from column (not working)
//    vec3 column_interaction_force = vec3(0,0,0);
//    float terrain_z = field_height_at(sphere.center.x, sphere.center.y);
//    if (min_dist < 1.5*smoke_layers[closest_layer_id].r)
//    {
//        vec3 force_vector = normalize(sphere.center - smoke_layers[closest_layer_id].center);
//        force_vector = vec3(force_vector.x, force_vector.y, 0);
//        column_interaction_force = 0.1 * m * (2*smoke_layers[closest_layer_id].r - min_dist) * force_vector;
//    }
//    else if (min_dist < 2.5*smoke_layers[closest_layer_id].r && terrain_z + 100. < sphere.center.z)
//    {
//        // horizontal friction
//        column_interaction_force = -0.001 * m/dt * vec3(sphere.speed.x, sphere.speed.y, 0);
//    }
//    if (smoke_layers[closest_layer_id].secondary_plume) column_interaction_force = vec3(0,0,0);
//    column_interaction_force = vec3(0,0,0);

    Eigen::Vector3f forces = gravity + buoyancy + friction;
    Eigen::Vector3f a = forces/m;
    Eigen::Vector3f v = sphere.speed + a * frame_nb*dt;
    Eigen::Vector3f p = sphere.center + v * frame_nb*dt;

    sphere.speed = v;
    sphere.center = p;

    if (!sphere.falling_disappeared) sphere_ground_collision(sphere, idx, frame_nb);
}

void scene_model::secondary_columns_creation()
{
    // merge close buffers
    for (unsigned int i = 0; i<falling_spheres_buffers.size(); i++)
    {
        for (unsigned int j = i+1; j<falling_spheres_buffers.size(); j++)
        {
            if ((falling_spheres_buffers[i][0].center - falling_spheres_buffers[j][0].center).norm() < 3 * falling_spheres_buffers[i][0].r)
            {
                for (unsigned int k = 0; k<falling_spheres_buffers[j].size(); k++)
                {
                    falling_spheres_buffers[i].push_back(falling_spheres_buffers[j][k]);
                }
                falling_spheres_buffers.erase(falling_spheres_buffers.begin()+j);
            }
        }
    }

    // find big enough buffers
    for (unsigned int i = 0; i<falling_spheres_buffers.size(); i++)
    {
        float wanted_ray = falling_spheres_buffers[i][0].r;
        float wanted_volume = wanted_ray * 3.14 * wanted_ray*wanted_ray;
        float sphere_volume = 4./3.*3.14*falling_spheres_buffers[i][0].r*falling_spheres_buffers[i][0].r*falling_spheres_buffers[i][0].r;
        float nb_spheres_needed = wanted_volume/sphere_volume;
        nb_spheres_needed = 6;

        //find closest layer
        Eigen::Vector3f center_i = falling_spheres_buffers[i][0].center;
        int closest_layer_id = smoke_layers.size()-1;
        float min_dist = (center_i - smoke_layers[closest_layer_id].center).norm();
        for (unsigned int j = 0; j<smoke_layers.size(); j++)
        {
            float dist = (center_i - smoke_layers[j].center).norm();
            if (dist < min_dist)
            {
                min_dist = dist;
                closest_layer_id = j;
            }
        }

        if (falling_spheres_buffers[i].size() > nb_spheres_needed*3 && min_dist > wanted_ray)
        {
            // emit layer
            add_smoke_layer(5, falling_spheres_buffers[i][0].rho, falling_spheres_buffers[i][0].charge, wanted_ray*0.75, falling_spheres_buffers[i][0].center, true);
            add_free_spheres_for_one_layer(smoke_layers.size()-1);
            total_layers_ejected++;
            //if (debug_mode) std::cout << "SECONDARY LAYER ADDED OK" << std::endl;

            // remove corresponding particles
            for (unsigned int j = 0; j<nb_spheres_needed; j++)
            {
                falling_spheres.push_back(falling_spheres_buffers[i][0]);
                falling_spheres[falling_spheres.size()-1].falling_disappeared = true;
                falling_spheres[falling_spheres.size()-1].rho = 10.;
                falling_spheres_buffers[i].erase(falling_spheres_buffers[i].begin());
            }
        }
    }

    // test to check closest spheres every timestep without buffers: too slow
//    for (unsigned int i = 0; i<falling_spheres.size(); i++)
//    {
//        if (falling_spheres[i].falling_under_atm_rho)
//        {
//            std::vector<unsigned int> close_particle_light;
//            for (unsigned int j = 0; j<falling_spheres.size(); j++)
//            {
//                if (falling_spheres[j].falling_under_atm_rho && norm(falling_spheres[i].center - falling_spheres[j].center)< 2*falling_spheres[i].r)
//                {
//                    close_particle_light.push_back(j);
//                }
//            }

//            vec3 center_i = falling_spheres[i].center;
//            //find closest layer
//            int closest_layer_id = smoke_layers.size()-1;
//            float min_dist = norm(center_i - smoke_layers[closest_layer_id].center);
//            for (unsigned int j = 0; j<smoke_layers.size(); j++)
//            {
//                float dist = norm (center_i - smoke_layers[j].center);
//                if (dist < min_dist)
//                {
//                    min_dist = dist;
//                    closest_layer_id = j;
//                }
//            }

//            if (close_particle_light.size() > 6 && min_dist > falling_spheres[i].r)
//            {
//                // emit layer
//                add_smoke_layer(5, falling_spheres[i].rho, falling_spheres[i].charge, falling_spheres[i].r, falling_spheres[i].center, true);
//                add_free_spheres_for_one_layer(smoke_layers.size()-1);
//                total_layers_ejected++;

//                // remove corresponding particles
//                for (unsigned int j = 0; j<6; j++)
//                {
//                    falling_spheres.erase(falling_spheres.begin() + close_particle_light[j]);
//                }
//            }
//        }
//    }

}

void scene_model::falling_spheres_update(unsigned int frame_nb)
{
    // edit spheres, attached or not to a buffer
    for (int i = falling_spheres.size()-1; i>=0; i--)
    {
        ground_falling_sphere_update(falling_spheres[i], i, frame_nb);
    }
    for (unsigned int i = 0; i<falling_spheres_buffers.size(); i++)
    {
        for (unsigned int j = 0; j<falling_spheres_buffers[i].size(); j++)
        {
            ground_falling_sphere_update(falling_spheres_buffers[i][j], -1, frame_nb);
        }
    }

    for (int i = falling_spheres.size()-1; i>=0; i--)
    {
        if (falling_spheres[i].falling_disappeared && falling_spheres[i].center.z() < -2000.)
        {
            falling_spheres.erase(falling_spheres.begin()+i);
        }
    }

    // emit new layers from buffers
    secondary_columns_creation();
}


//------------------------------------------------------------
//---------------------- FREE SPHERES ------------------------
//------------------------------------------------------------ */

void scene_model::add_free_sphere(unsigned int i, float angle, float size_fac)
{
    free_sphere_params sphere(smoke_layers[i].center, angle, size_fac * smoke_layers[i].r, smoke_layers[i].v.z(), smoke_layers[i].charge);
    sphere.size_factor = size_fac;
    sphere.rho = smoke_layers[i].rho;
    float volume = 4.0/3.0 * 3.14 * sphere.r*sphere.r*sphere.r;
    sphere.mass = sphere.rho/volume;
    sphere.closest_layer_idx = i;

    if (smoke_layers[i].secondary_plume)
    {
        sphere.secondary_column = true;
        sphere.closest_layer_idx = i;
    }

    for (unsigned int i = 0; i<subspheres_number; i++)
    {
        // random angles
        float theta = 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float phi = 2*3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        subsphere_params subs = subsphere_params();
        subs.parent_id = free_spheres.size();
        subs.relative_position = Eigen::Vector3f(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
        subs.center = sphere.center + sphere.r * subs.relative_position;
        //subs.size_ratio = 0.10 + 0.25 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        subs.size_ratio = 0.2;
        subs.r = sphere.r*subs.size_ratio;
        subs.charge = sphere.charge;
        subs.speed = sphere.speed;

        for (unsigned int j = 0; j<subsubspheres_number; j++)
        {
            float theta2 = 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            float phi2 = 2*3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

            subsphere_params subsubs = subsphere_params();
            subsubs.parent_id = s2_spheres.size();
            subsubs.relative_position = Eigen::Vector3f(sin(theta2)*cos(phi2), sin(theta2)*sin(phi2), cos(theta2));
            subsubs.center = subs.center + subs.r * subsubs.relative_position;
            subsubs.r = subs.r/5.;
            subsubs.charge = subs.charge;
            subsubs.speed = subs.speed;
            s3_spheres.push_back(subsubs);
        }
        s2_spheres.push_back(subs);
    }
    free_spheres.push_back(sphere);
}

void scene_model::add_free_spheres_for_one_layer(unsigned int i)
{
    unsigned int nb_spheres = 6;
    float angle_offset = 2*3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

    //determine size differences
    std::vector<float> sizes;
    for (unsigned int k = 0; k<nb_spheres; k++)
    {
        sizes.push_back(0.5f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
    }
    float total = sizes[0] + sizes[1] + sizes[2] + sizes[3] + sizes[4] + sizes[5];

    float factor = (float)nb_spheres / (1.0f * total); // JG PARAM - global scaling factor to make spheres smaller

    for (unsigned int k = 0; k<sizes.size(); k++)
    {
        sizes[k] *= factor;
    }

    // add spheres
    for (unsigned int j = 0; j<nb_spheres; j++)
    {
        float angle = (float)j*2.0*3.14/6.0;
        add_free_sphere(i, angle+angle_offset, sizes[j]);
    }
}

float compute_gaussian_speed_in_layer(float v_z, float max_r, float r)
{
    float A = max_r/(2.*sqrt(log(2)));
    return 2. * v_z * exp(-r*r/(A*A));
}

void scene_model::subdivide_and_make_falling(unsigned int i)
{
    // update subspheres
    update_subspheres_params();

    // define ray of falling spheres; density same as free sphere
    float falling_ray = free_spheres[i].r/5.;
    float falling_volume = 4./3.*3.14*falling_ray*falling_ray*falling_ray;
    float free_sphere_volume = 4./3.*3.14*free_spheres[i].r*free_spheres[i].r*free_spheres[i].r;
    float n_float = free_sphere_volume/falling_volume;
    int n = (int)n_float;

    // start adding subspheres as falling spheres
    for (unsigned int j = 0; j<s2_spheres.size(); j++)
    {
        if (s2_spheres[j].parent_id == i)
        {
            free_sphere_params sphere = free_sphere_params(s2_spheres[j].center, s2_spheres[j].r, free_spheres[i].rho, free_spheres[i].charge);
            sphere.falling = true;
            sphere.stagnate = false;
            falling_spheres.push_back(sphere);
            if (n>=0) n--;
        }
    }

    // add all falling spheres
    if (n>0)
    {
        for (unsigned int j = 0; j<n; j++)
        {
            float rand_r = free_spheres[i].r * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            float rand_phi = 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            float rand_theta = 2*3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            Eigen::Vector3f rand_vec(sin(rand_theta)*cos(rand_phi), sin(rand_theta)*sin(rand_phi), cos(rand_theta));
            Eigen::Vector3f new_center = free_spheres[i].center + rand_r * rand_vec; //random position inside free sphere
            free_sphere_params sphere = free_sphere_params(new_center, falling_ray, free_spheres[i].rho, free_spheres[i].charge);
            falling_spheres.push_back(sphere);
        }
    }

    // make free sphere falling
    free_spheres[i].falling = true;
    if (debug_mode) std::cout << i << " falls" << std::endl;
}

void scene_model::update_free_spheres()
{
    int size = free_spheres.size() - 1;
//#pragma omp parallel for
    for (int i = size; i>=0; --i)
    {
        free_sphere_params& sphere_i = free_spheres[i];

        if (!sphere_i.stagnate_long && !sphere_i.falling)
        {
            // identify closest layer which is rising or begins falling (later : all layers in which the sphere is)
            int closest_layer_id = smoke_layers.size()-1;
            float min_dist = (sphere_i.center - smoke_layers[closest_layer_id].center).norm();
            for (unsigned int j = 0; j<smoke_layers.size(); j++)
            {
                //float dist = abs(free_spheres[i].center.z - smoke_layers[j].center.z);
                float dist = (sphere_i.center - smoke_layers[j].center).norm();
                if (dist < min_dist && (smoke_layers[j].rising || smoke_layers[j].begin_falling) && !smoke_layers[j].stagnates_long)
                {
                    min_dist = dist;
                    closest_layer_id = j;
                }
            }
            if (sphere_i.secondary_column) closest_layer_id = sphere_i.closest_layer_idx;
            if (sphere_i.stagnate || sphere_i.stagnate_long) closest_layer_id = sphere_i.closest_layer_idx;
            closest_layer_id = sphere_i.closest_layer_idx;

            // check if closest layer begins falling (if so, make sphere falling)
            if (smoke_layers[closest_layer_id].begin_falling)
            {
                subdivide_and_make_falling(i);
            }
            else if (smoke_layers[closest_layer_id].stagnates && !smoke_layers[closest_layer_id].stagnates_long && !sphere_i.stagnate)
            {
                // check if densities have become equal: if so, keep altitude in memory
                sphere_i.stagnate = true;
                sphere_i.closest_layer_idx = closest_layer_id;
                sphere_i.stagnation_altitude = sphere_i.center.z();
            }
            else if (sphere_i.stagnate && !sphere_i.stagnate_long && smoke_layers[sphere_i.closest_layer_idx].stagnates_long)
            {
                // check if closest layer has reached max altitude: if so, keep sphere position in memory for stagnation spreading function
                sphere_i.stagnate_long = true;
                sphere_i.max_altitude = sphere_i.center.z();
                sphere_i.center_at_max_altitude = sphere_i.center;
                sphere_i.layer_center_at_max_altitude = smoke_layers[closest_layer_id].center;
                sphere_i.xy_at_max_altitude = sqrt(sphere_i.center.x() * sphere_i.center.x() + sphere_i.center.y() * sphere_i.center.y());
            }
            else
            {
                // get axial and radial composants relative to layer center
                Eigen::Vector3f p_relative = sphere_i.center - smoke_layers[closest_layer_id].center;
                Eigen::Vector3f p_axial = p_relative.dot(smoke_layers[closest_layer_id].plume_axis) * smoke_layers[closest_layer_id].plume_axis;
                Eigen::Vector3f p_radial = p_relative - p_axial;

                // update rotation axis with new radial vector (can change upon time because axis changes)
                sphere_i.angle_vector = p_radial/(p_radial).norm();
                //sphere_i.rotation_axis = normalize(cross(smoke_layers[closest_layer_id].plume_axis, p_radial));

                // update radial position if too close from plume axis
                //if (norm(p_radial) < smoke_layers[closest_layer_id].r)*/ free_spheres[i].center = smoke_layers[closest_layer_id].center + p_axial + smoke_layers[closest_layer_id].r * normalize(p_radial);

                // update size according to layer
                // size of layer + perturbation (some spheres should grow much more than others, to create diversity)
                float new_r = sphere_i.size_factor * smoke_layers[closest_layer_id].r;

                // update radial speed and position by adding perturbation (one part is random and one depends on radial position, so that spheres do not stay in the middle of the plume and do not go away)
                float random_f = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                sphere_i.perturbation += new_r * 0.001 * 2.0 * (random_f - 0.5);
                if (sphere_i.perturbation > 100.) sphere_i.perturbation = 100.;
                if (sphere_i.perturbation < -100.) sphere_i.perturbation = -100.;
                float new_speed = compute_gaussian_speed_in_layer(smoke_layers[closest_layer_id].speed_along_axis, 2.*smoke_layers[closest_layer_id].r, (p_relative).norm()); // axial speed
                if (smoke_layers[closest_layer_id].stagnates) new_speed = smoke_layers[closest_layer_id].speed_along_axis;
                float radial_speed = (smoke_layers[closest_layer_id].r - (p_radial).norm())/0.5;
                if (smoke_layers[closest_layer_id].stagnates) radial_speed = 0;
                sphere_i.perturbation += radial_speed;
                //if (smoke_layers[closest_layer_id].stagnates) std::cout << sphere_i.perturbation << std::endl;
                if (smoke_layers[closest_layer_id].stagnates && sphere_i.perturbation < 0) sphere_i.perturbation = 0;
                radial_speed = 0;
                //std::cout << norm(p_relative) - smoke_layers[closest_layer_id].r << " " << free_spheres[i].perturbation << " " << radial_speed << std::endl;

                // update rotation speed with new speed and ray
                float new_angular_speed = new_speed/new_r;

                // update
                sphere_i.speed = new_speed * smoke_layers[closest_layer_id].plume_axis + (sphere_i.perturbation) * sphere_i.angle_vector;
                sphere_i.r = new_r;
                sphere_i.relative_distance = (sphere_i.center-smoke_layers[closest_layer_id].center).norm();
                //sphere_i.density_loss_last_dt += sphere_i.rho - smoke_layers[closest_layer_id].rho;
                sphere_i.rho = smoke_layers[closest_layer_id].rho;

                if (!sphere_i.stagnate && smoke_layers[closest_layer_id].theta >1)
                {
                    sphere_i.angular_speed = new_angular_speed;
                    float new_angle = sphere_i.current_angle + sphere_i.angular_speed * dt;
                    sphere_i.current_angle = new_angle;
                }
                else sphere_i.angular_speed = 0;
            }
        }
    }
//#pragma omp barrier
    
#pragma omp parallel for
    // update positions
    for (unsigned int i = 0; i<free_spheres.size(); ++i)
    {
        if (!free_spheres[i].stagnate_long && !free_spheres[i].falling)
        {
            free_spheres[i].center += dt * free_spheres[i].speed;
        }
        else if (free_spheres[i].falling)
        {
            free_spheres[i].center.z() = smoke_layers[free_spheres[i].closest_layer_idx].center.z();
        }
        //update_spheres_on_free_sphere(i);
    }
#pragma omp barrier
    
    // make begin_falling layers falling
    for (int i = smoke_layers.size()-1; i>=0; i--)
    {
        if (smoke_layers[i].begin_falling && smoke_layers.size() >= 1)
        {
            smoke_layers[i].falling = true;
            smoke_layers[i].begin_falling = false;
            smoke_layers[i].rising = false;
            smoke_layers[i].plume = false;
            //smoke_layers[i].center.z = -2000;
        }
    }
}


//------------------------------------------------------------
//---------------------- STAGNATION --------------------------
//------------------------------------------------------------ */

void scene_model::update_stagnation_spheres()
{
    int size = free_spheres.size() - 1;
    //#pragma omp parallel for
    for (int i = size; i>=0; --i)
    {
        if (free_spheres[i].stagnate_long)
        {
            // closest layer
            unsigned int closest_layer_id = free_spheres[i].closest_layer_idx;

            // get wind
            //Eigen::Vector3f wind = compute_wind_vector(free_spheres[i].center.z());
            Eigen::Vector2f wind_2d = scene_model::sampleVel(free_spheres[i].center);
            Eigen::Vector3f wind(wind_2d.x(), wind_2d.y(), 0.f);

            // get radial composant relative to layer center
            Eigen::Vector3f p_radial = Eigen::Vector3f(free_spheres[i].center.x(), free_spheres[i].center.y(), 0);

            // update rotation axis with new radial vector (can change upon time because axis changes)
            free_spheres[i].angle_vector = p_radial.normalized();

            // update radial speed and position by adding perturbation
            float speed_factor = (Eigen::Vector3f(free_spheres[i].center.x(), free_spheres[i].center.y(), 0)).norm() / 2000.;
            free_spheres[i].perturbation = stagnation_speed/speed_factor;
            if (wind.norm() != 0)
            {
                speed_factor = p_radial.norm() / 2000.;
                free_spheres[i].perturbation = stagnation_speed/speed_factor;
                if (p_radial.dot(free_spheres[i].speed) < 0) free_spheres[i].perturbation = 0;
            }

            //free_spheres[i].angle_vector = normalize(vec3(free_spheres[i].angle_vector.x, free_spheres[i].angle_vector.y, 0));

            // update
            free_spheres[i].speed = 0*free_spheres[i].speed + 1*(free_spheres[i].perturbation) * free_spheres[i].angle_vector;

            float dxy = (free_spheres[i].speed).norm() * dt;
            Eigen::Vector3f relative_pos_at_max_alt = free_spheres[i].center_at_max_altitude - free_spheres[i].layer_center_at_max_altitude;
            float relative_xy_at_max_alt = sqrt(relative_pos_at_max_alt.x() * relative_pos_at_max_alt.x() + relative_pos_at_max_alt.y() * relative_pos_at_max_alt.y());
            float a = (free_spheres[i].max_altitude - free_spheres[i].stagnation_altitude) * relative_xy_at_max_alt;

            Eigen::Vector3f relative_pos = free_spheres[i].center - free_spheres[i].layer_center_at_max_altitude;
            float relative_xy = sqrt(relative_pos.x() * relative_pos.x() + relative_pos.y() * relative_pos.y());
            float dz = -a * dxy / ((relative_xy)*(relative_xy));
            if (relative_xy < relative_xy_at_max_alt) dz = 0;
            float new_z = free_spheres[i].stagnation_altitude + a / (relative_xy);
            if (relative_xy <= relative_xy_at_max_alt) new_z = free_spheres[i].max_altitude;

            free_spheres[i].center.z() = new_z;
            free_spheres[i].angular_speed = 0;
            free_spheres[i].density_loss_last_dt += free_spheres[i].rho - smoke_layers[closest_layer_id].rho;
            free_spheres[i].rho = smoke_layers[closest_layer_id].rho;

            if (free_spheres[i].z_factor > 0.5 /*- 0.003*norm(wind)*/) free_spheres[i].z_factor -= 0.005* dt + 0.0001 * wind.norm() * dt;

            free_spheres[i].speed = wind;
        }
    }

    // update positions
#pragma omp parallel for
    for (unsigned int i = 0; i<free_spheres.size(); i++)
    {
        if (free_spheres[i].stagnate_long)
        {
            free_spheres[i].center += dt * free_spheres[i].speed;
            //update_spheres_on_stagnation_sphere(i);
            if (free_spheres[i].closest_layer_idx != -1)
            {
                Eigen::Vector3f relative_position = free_spheres[i].center - smoke_layers[free_spheres[i].closest_layer_idx].center;
                free_spheres[i].relative_distance = relative_position.norm();
                free_spheres[i].angle_vector = relative_position.normalized();
            }
        }
    }
#pragma omp barrier
}


//------------------------------------------------------------
//------------------------- EXPORT ---------------------------
//------------------------------------------------------------ */

void scene_model::update_subspheres_params()
{
    // normal spheres
    for (size_t i = 0; i<s2_spheres.size(); i++)
    {
        // apply parent transfo : s_i = s/5, t_i = s*R*offset_i  + t
        unsigned int idx_parent = s2_spheres[i].parent_id;

        if (!free_spheres[idx_parent].falling)
        {
            Eigen::Matrix3f const R_parent = (Eigen::AngleAxisf(free_spheres[idx_parent].current_angle, free_spheres[idx_parent].rotation_axis)).toRotationMatrix();
            float r_parent = free_spheres[idx_parent].r;
            Eigen::Vector3f t_parent = free_spheres[idx_parent].center;

            s2_spheres[i].r = r_parent*s2_spheres[i].size_ratio;
            s2_spheres[i].center = r_parent*R_parent*s2_spheres[i].relative_position + t_parent;

            if (free_spheres[idx_parent].stagnate_long)
            {
                float z_distance_from_center = s2_spheres[i].center.z() - t_parent.z();
                s2_spheres[i].z_factor = free_spheres[idx_parent].z_factor;
                s2_spheres[i].center.z() = t_parent.z() + s2_spheres[i].z_factor * z_distance_from_center;
            }

            Eigen::Vector3f tangent = free_spheres[idx_parent].rotation_axis.cross(s2_spheres[i].center - t_parent);
            tangent.normalize();
            s2_spheres[i].speed = free_spheres[idx_parent].speed + r_parent*free_spheres[idx_parent].angular_speed*tangent;
        }
    }

    for (size_t j = 0; j<s3_spheres.size(); j++)
    {
        // apply parent transfo : sj = si/5, tj = si*offsetj  + ti
        unsigned int idx_parent = s3_spheres[j].parent_id;
        if (!free_spheres[s2_spheres[idx_parent].parent_id].falling)
        {
            float r_parent = s2_spheres[idx_parent].r;
            Eigen::Vector3f t_parent = s2_spheres[idx_parent].center;

            s3_spheres[j].r = r_parent/5.;
            s3_spheres[j].center = r_parent*s3_spheres[j].relative_position + t_parent;
        }
    }
}

/*
void scene_model::export_spheres_combined(unsigned int plume_export_frequency)
{
    update_subspheres_params();

    unsigned int frame_nb = frame_count/plume_export_frequency;

    std::string spheres_centers_file;
    std::string densities_file;
   
    spheres_centers_file = "../output/plumespheres_" + std::to_string(frame_nb) + ".bin";
    std::ofstream spheres_centers = std::ofstream(spheres_centers_file.c_str());

    if (spheres_centers)
    {
        // total sphere count
        spheres_centers << (int) (free_spheres.size() + s2_spheres.size() + falling_spheres.size() + falling_spheres_buffers.size()) << std::endl;
    }
    
    // free + stagnation spheres NEW
    for (unsigned int i = 0; i<free_spheres.size(); i++)
    {
        free_sphere_params s = free_spheres[i];
        if (spheres_centers)
            spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << " " << s.rho << std::endl;
    }

    // subspheres NEW
    for (unsigned int i = 0; i<s2_spheres.size(); i++)
    {
        unsigned int parent_idx = s2_spheres[i].parent_id;
        if (!free_spheres[parent_idx].falling)
        {
            subsphere_params s = s2_spheres[i];
            if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << " " << free_spheres[parent_idx].rho << std::endl;
        }
    }
   
    // falling spheres
    for (unsigned int i = 0; i<falling_spheres.size(); i++)
    {
        free_sphere_params s = falling_spheres[i];
        if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << " " << s.rho << std::endl;
    }

    // falling spheres buffers
    for (unsigned int k = 0; k<falling_spheres_buffers.size(); k++)
    {
        for (unsigned int i = 0; i<falling_spheres_buffers[k].size(); i++)
        {
            free_sphere_params s = falling_spheres_buffers[k][i];
            if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << " " << s.rho << std::endl;
        }
    }

}*/
            
void scene_model::export_ashfall(unsigned int plume_export_frequency)
{
    //update_subspheres_params();

    unsigned int frame_nb = frame_count/plume_export_frequency;

    std::string spheres_x_file;
    std::string spheres_y_file;
    std::string spheres_z_file;
    std::string spheres_r_file;
    std::string ash_loss_file;
    if (frame_nb < 10)
    {
        spheres_x_file = "../output/spheres_x_00" + std::to_string(frame_nb) + ".txt";
        spheres_y_file = "../output/spheres_y_00" + std::to_string(frame_nb) + ".txt";
        spheres_z_file = "../output/spheres_z_00" + std::to_string(frame_nb) + ".txt";
        spheres_r_file = "../output/spheres_r_00" + std::to_string(frame_nb) + ".txt";
        ash_loss_file = "../output/ash_loss_00" + std::to_string(frame_nb) + ".txt";
    }
    else if (frame_nb < 100)
    {
        spheres_x_file = "../output/spheres_x_0" + std::to_string(frame_nb) + ".txt";
        spheres_y_file = "../output/spheres_y_0" + std::to_string(frame_nb) + ".txt";
        spheres_z_file = "../output/spheres_z_0" + std::to_string(frame_nb) + ".txt";
        spheres_r_file = "../output/spheres_r_0" + std::to_string(frame_nb) + ".txt";
        ash_loss_file = "../output/ash_loss_0" + std::to_string(frame_nb) + ".txt";
    }
    else
    {
        spheres_x_file = "../output/spheres_x_" + std::to_string(frame_nb) + ".txt";
        spheres_y_file = "../output/spheres_y_" + std::to_string(frame_nb) + ".txt";
        spheres_z_file = "../output/spheres_z_" + std::to_string(frame_nb) + ".txt";
        spheres_r_file = "../output/spheres_r_" + std::to_string(frame_nb) + ".txt";
        ash_loss_file = "../output/ash_loss_" + std::to_string(frame_nb) + ".txt";
    }
    std::ofstream spheres_x = std::ofstream(spheres_x_file.c_str());
    std::ofstream spheres_y = std::ofstream(spheres_y_file.c_str());
    std::ofstream spheres_z = std::ofstream(spheres_z_file.c_str());
    std::ofstream spheres_r = std::ofstream(spheres_r_file.c_str());
    std::ofstream ash_loss = std::ofstream(ash_loss_file.c_str());


    // free + stagnation spheres NEW
    for (unsigned int i = 0; i<free_spheres.size(); i++)
    {
        if (true)
        {
            free_sphere_params s = free_spheres[i];
            if (spheres_x) spheres_x << s.center.x()<< std::endl;
            if (spheres_y) spheres_y << s.center.y()<< std::endl;
            if (spheres_z) spheres_z << s.center.z()<< std::endl;
            if (spheres_r) spheres_r << s.r << std::endl;
            if (ash_loss) ash_loss << s.density_loss_last_dt << std::endl;
        }
    }

}

void scene_model::export_spheres(unsigned int plume_export_frequency)
{
    update_subspheres_params();

    unsigned int frame_nb = frame_count/plume_export_frequency;

    //std::string s1_spheres_file;
    //std::string s1s2_spheres_file;

    std::string spheres_centers_file;
    std::string densities_file;
    if (frame_nb < 10)
    {
        spheres_centers_file = "../output/spheres_00" + std::to_string(frame_nb) + ".txt";
        densities_file = "../output/densities_00" + std::to_string(frame_nb) + ".txt";
        //s1_spheres_file = "../output/s1_spheres_00" + std::to_string(frame_nb) + ".txt";
        //s1s2_spheres_file = "../output/s1s2_spheres_00" + std::to_string(frame_nb) + ".txt";
    }
    else if (frame_nb < 100)
    {
        spheres_centers_file = "../output/spheres_0" + std::to_string(frame_nb) + ".txt";
        densities_file = "../output/densities_0" + std::to_string(frame_nb) + ".txt";
        //s1_spheres_file = "../output/s1_spheres_0" + std::to_string(frame_nb) + ".txt";
        //s1s2_spheres_file = "../output/s1s2_spheres_0" + std::to_string(frame_nb) + ".txt";
    }
    else
    {
        spheres_centers_file = "../output/spheres_" + std::to_string(frame_nb) + ".txt";
        densities_file = "../output/densities_" + std::to_string(frame_nb) + ".txt";
        //s1_spheres_file = "../output/s1_spheres_" + std::to_string(frame_nb) + ".txt";
        //s1s2_spheres_file = "../output/s1s2_spheres_" + std::to_string(frame_nb) + ".txt";
    }
    std::ofstream spheres_centers = std::ofstream(spheres_centers_file.c_str());
    std::ofstream densities = std::ofstream(densities_file.c_str());

    //std::ofstream s1_spheres = std::ofstream(s1_spheres_file.c_str());
    //std::ofstream s1s2_spheres = std::ofstream(s1s2_spheres_file.c_str());

    // free + stagnation spheres NEW
    for (unsigned int i = 0; i<free_spheres.size(); i++)
    {
        if (true)
        {
            free_sphere_params s = free_spheres[i];
            //if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << std::endl;
            //if (densities) densities << s.rho << std::endl;

            //if (s1_spheres) s1_spheres << s.center.x/1000. << " " << s.center.z/1000. << " " << s.center.y/1000. << " " << s.r/1000. << std::endl;
            //if (s1s2_spheres) s1s2_spheres << s.center.x/1000. << " " << s.center.z/1000. << " " << s.center.y/1000. << " " << s.r/1000. << std::endl;
        }
    }

    // subspheres NEW
    for (unsigned int i = 0; i<s2_spheres.size(); i++)
    {
        unsigned int parent_idx = s2_spheres[i].parent_id;
        if (!free_spheres[parent_idx].falling)
        {
            subsphere_params s = s2_spheres[i];
            if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << std::endl;
            if (densities) densities << free_spheres[parent_idx].rho << std::endl;
            //if (s1s2_spheres) s1s2_spheres << s.center.x/1000. << " " << s.center.z/1000. << " " << s.center.y/1000. << " " << s.r/1000. << std::endl;
        }
    }
    // uncomment to get 3rd LoD (takes much longer + files become much heavier)
//    for (unsigned int i = 0; i<s3_spheres.size(); i++)
//    {
//        unsigned int parent_idx = s3_spheres[i].parent_id;
//        if (!free_spheres[s2_spheres[parent_idx].parent_id].falling)
//        {
//            subsphere_params s = s3_spheres[i];
//            if (spheres_centers) spheres_centers << s.center.x/1000. << " " << s.center.z/1000. << " " << s.center.y/1000. << " " << s.r/1000. << std::endl;
//            if (densities) densities << free_spheres[s2_spheres[parent_idx].parent_id].rho << std::endl;
//        }
//    }

    // falling spheres
    for (unsigned int i = 0; i<falling_spheres.size(); i++)
    {
        free_sphere_params s = falling_spheres[i];
        if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << std::endl;
        if (densities) densities << s.rho << std::endl;

        //spheres on spheres
//        for (unsigned int j = 0; j<s.spheres_on_sphere.size(); j++)
//        {
//            sphere_on_free_spheres_params s_2 = s.spheres_on_sphere[j];
//            if (spheres_centers) spheres_centers << s_2.center.x/1000. << " " << s_2.center.z/1000. << " " << s_2.center.y/1000. << " " << s_2.r/1000. << std::endl;
//            if (densities) densities << s.rho << std::endl;

//            //spheres on spheres on spheres
//            for (unsigned int k = 0; k<s_2.spheres_on_sphere.size(); k++)
//            {
//                sphere_on_free_spheres_params s_3 = s_2.spheres_on_sphere[k];
//                if (spheres_centers) spheres_centers << s_3.center.x/1000. << " " << s_3.center.z/1000. << " " << s_3.center.y/1000. << " " << s_3.r/1000. << std::endl;
//                if (densities) densities << s.rho << std::endl;


//            }
//        }
    }

    // falling spheres buffers
    for (unsigned int k = 0; k<falling_spheres_buffers.size(); k++)
    {
        for (unsigned int i = 0; i<falling_spheres_buffers[k].size(); i++)
        {
            free_sphere_params s = falling_spheres_buffers[k][i];
            if (spheres_centers) spheres_centers << s.center.x()/1000. << " " << s.center.z()/1000. << " " << s.center.y()/1000. << " " << s.r/1000. << std::endl;
            if (densities) densities << s.rho << std::endl;

            //spheres on spheres
//            for (unsigned int j = 0; j<s.spheres_on_sphere.size(); j++)
//            {
//                sphere_on_free_spheres_params s_2 = s.spheres_on_sphere[j];
//                if (spheres_centers) spheres_centers << s_2.center.x/1000. << " " << s_2.center.z/1000. << " " << s_2.center.y/1000. << " " << s_2.r/1000. << std::endl;
//                if (densities) densities << s.rho << std::endl;

//                //spheres on spheres on spheres
//                for (unsigned int k = 0; k<s_2.spheres_on_sphere.size(); k++)
//                {
//                    sphere_on_free_spheres_params s_3 = s_2.spheres_on_sphere[k];
//                    if (spheres_centers) spheres_centers << s_3.center.x/1000. << " " << s_3.center.z/1000. << " " << s_3.center.y/1000. << " " << s_3.r/1000. << std::endl;
//                    if (densities) densities << s.rho << std::endl;


//                }
//            }
        }
    }

}

void scene_model::export_spheres_binary(unsigned int plume_export_frequency, bool headercount)
{
    update_subspheres_params();

    unsigned int frame_nb = frame_count/plume_export_frequency;
    float minc = 100000.0f, maxc = 0.0f, meanc = 0.0f;
    int countc = 0;

    std::string spheres_centers_file;
    std::string densities_file;
    std::string velocity_file;
    std::string charge_file;
    std::string type_file;
    
    if (frame_nb < 10)
    {
        spheres_centers_file = "../output/spheres_00" + std::to_string(frame_nb) + ".dat";
        densities_file = "../output/densities_00" + std::to_string(frame_nb) + ".dat";
        velocity_file = "../output/velocity_00" + std::to_string(frame_nb) + ".dat";
        charge_file = "../output/charge_00" + std::to_string(frame_nb) + ".dat";
        type_file = "../output/type_00" + std::to_string(frame_nb) + ".dat";
    }
    else if (frame_nb < 100)
    {
        spheres_centers_file = "../output/spheres_0" + std::to_string(frame_nb) + ".dat";
        densities_file = "../output/densities_0" + std::to_string(frame_nb) + ".dat";
        velocity_file = "../output/velocity_0" + std::to_string(frame_nb) + ".dat";
        charge_file = "../output/charge_0" + std::to_string(frame_nb) + ".dat";
        type_file = "../output/type_0" + std::to_string(frame_nb) + ".dat";
    }
    else
    {
        spheres_centers_file = "../output/spheres_" + std::to_string(frame_nb) + ".dat";
        densities_file = "../output/densities_" + std::to_string(frame_nb) + ".dat";
        velocity_file = "../output/velocity_" + std::to_string(frame_nb) + ".dat";
        charge_file = "../output/charge_" + std::to_string(frame_nb) + ".dat";
        type_file = "../output/type_" + std::to_string(frame_nb) + ".dat";
    }
    std::ofstream spheres_centers(spheres_centers_file, std::ofstream::binary | std::ofstream::out);
    std::ofstream densities(densities_file, std::ofstream::binary | std::ofstream::out);
    std::ofstream velocity(velocity_file, std::ofstream::binary | std::ofstream::out);
    std::ofstream charge(charge_file, std::ofstream::binary | std::ofstream::out);
    std::ofstream sphtype(type_file, std::ofstream::binary | std::ofstream::out);


    // out.write(reinterpret_cast<const char*>( &f ), sizeof( float ));
    float x=0, y=0, z=0, r=0;
    float vx=0, vy=0, vz=0, chg=0;
    
    if(headercount && spheres_centers)
    {
        int cnt = (int) (/*free_spheres.size() +*/ s2_spheres.size() + falling_spheres.size() + falling_spheres_buffers.size());
        spheres_centers.write(reinterpret_cast<const char*>( & cnt ), sizeof( int ));
    }
        
    float minr = 100000.0f, maxr = 0.0f;
    float minx = 100000.0f, maxx = -100000.0f, miny = 100000.0f, maxy = -100000.0f, minz = 100000.0f, maxz = -100000.0f;
    
    // free + stagnation spheres
    /*for (unsigned int i = 0; i<free_spheres.size(); i++)
    {
        if (true)
        {
            free_sphere_params s = free_spheres[i];
            x = s.center.x()/1000.;
            y = s.center.y()/1000.;
            z = s.center.z()/1000.;
            r = s.r/1000.;
            
            if(s.center.x() > maxx)
                maxx = s.center.x();
            if(s.center.y() > maxy)
                maxy = s.center.y();
            if(s.center.z() > maxz)
                maxz = s.center.z();
            if(s.center.x() < minx)
                minx = s.center.x();
            if(s.center.y() < miny)
                miny = s.center.y();
            if(s.center.z() < minz)
                minz = s.center.z();
            if(s.r > maxr)
                maxr = s.r;
            if(s.r < minr)
                minr = s.r;
            
            if (spheres_centers)
            {
                spheres_centers.write(reinterpret_cast<const char*>( &x ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &z ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &y ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &r ), sizeof( float ));
            }
            //std::string center_str = std::to_string(s.center.x()/1000.) + " " + std::to_string(s.center.z()/1000.) + " " + std::to_string(s.center.y()/1000.) + " " + std::to_string(s.r/1000.);
            //const char* center_char = center_str.c_str();
            //if (spheres_centers) spheres_centers.write(center_char, sizeof(center_char));
            //const char* rho = std::to_string(s.rho).c_str();
            float rho = s.rho;
            if (densities) densities.write(reinterpret_cast<const char*>( &rho ), sizeof( float ));

            vx = s.speed.x();
            vy = s.speed.y();
            vz = s.speed.z();
            if (velocity)
            {
                velocity.write(reinterpret_cast<const char*>( &vx ), sizeof( float ));
                velocity.write(reinterpret_cast<const char*>( &vz ), sizeof( float ));
                velocity.write(reinterpret_cast<const char*>( &vy ), sizeof( float ));
            }
            
            chg = s.charge;
            if(charge)
                charge.write(reinterpret_cast<const char*>( &chg ), sizeof( float ));
            if(chg < minc)
                minc = chg;
            if(chg > maxc)
                maxc = chg;
            meanc += chg;
            countc++;
            
            int val = 0;
            if(sphtype)
                sphtype.write(reinterpret_cast<const char*>( &val ), sizeof( int ));
        }
    }*/
    /*
    std::cerr << "FREE AND STAGNATION SPHERES" << std::endl;
    std::cerr << "maxr = " << maxr << " minr = " << minr << std::endl;
    std::cerr << "maxx = " << maxx << " minx = " << minx << std::endl;
    std::cerr << "maxy = " << maxy << " miny = " << miny << std::endl;
    std::cerr << "maxz = " << maxz << " minz = " << minz << std::endl;
     */

    // subspheres
    minr = 100000.0f; maxr = 0.0f;
    minx = 100000.0f; maxx = -100000.0f; miny = 100000.0f; maxy = -100000.0f; minz = 100000.0f; maxz = -100000.0f;
    for (unsigned int i = 0; i<s2_spheres.size(); i++)
    {
        unsigned int parent_idx = s2_spheres[i].parent_id;
        if (!free_spheres[parent_idx].falling)
        {
            subsphere_params s = s2_spheres[i];
            x = s.center.x()/1000.;
            y = s.center.y()/1000.;
            z = s.center.z()/1000.;
            r = s.r/1000.;
            
            if(s.center.x() > maxx)
                maxx = s.center.x();
            if(s.center.y() > maxy)
                maxy = s.center.y();
            if(s.center.z() > maxz)
                maxz = s.center.z();
            if(s.center.x() < minx)
                minx = s.center.x();
            if(s.center.y() < miny)
                miny = s.center.y();
            if(s.center.z() < minz)
                minz = s.center.z();
            if(s.r > maxr)
                maxr = s.r;
            if(s.r < minr)
                minr = s.r;
            
            if (spheres_centers)
            {
                spheres_centers.write(reinterpret_cast<const char*>( &x ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &z ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &y ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &r ), sizeof( float ));
            }
            //std::string center_str = std::to_string(s.center.x()/1000.) + " " + std::to_string(s.center.z()/1000.) + " " + std::to_string(s.center.y()/1000.) + " " + std::to_string(s.r/1000.);
            //const char* center_char = center_str.c_str();
            //if (spheres_centers) spheres_centers.write(center_char, sizeof(center_char));
            //const char* rho = std::to_string(free_spheres[parent_idx].rho).c_str();
            float rho = free_spheres[parent_idx].rho;
            if (densities) densities.write(reinterpret_cast<const char*>( &rho ), sizeof( float ));

            vx = s.speed.x();
            vy = s.speed.y();
            vz = s.speed.z();
            if (velocity)
            {
                velocity.write(reinterpret_cast<const char*>( &vx ), sizeof( float ));
                velocity.write(reinterpret_cast<const char*>( &vz ), sizeof( float ));
                velocity.write(reinterpret_cast<const char*>( &vy ), sizeof( float ));
            }
            
            chg = s.charge;
            if(charge)
                charge.write(reinterpret_cast<const char*>( &chg ), sizeof( float ));
            
            int val = 1;
            if(sphtype)
                sphtype.write(reinterpret_cast<const char*>( &val ), sizeof( int ));
        }
    }
    
    /*
    std::cerr << "SUBSPHERES" << std::endl;
    std::cerr << "maxr = " << maxr << " minr = " << minr << std::endl;
    std::cerr << "maxx = " << maxx << " minx = " << minx << std::endl;
    std::cerr << "maxy = " << maxy << " miny = " << miny << std::endl;
    std::cerr << "maxz = " << maxz << " minz = " << minz << std::endl;
     */

    // falling spheres
    minr = 100000.0f; maxr = 0.0f;
    minx = 100000.0f; maxx = -100000.0f; miny = 100000.0f; maxy = -100000.0f; minz = 100000.0f; maxz = -100000.0f;
    for (unsigned int i = 0; i<falling_spheres.size(); i++)
    {
        free_sphere_params s = falling_spheres[i];
        x = s.center.x()/1000.;
        y = s.center.y()/1000.;
        z = s.center.z()/1000.;
        r = s.r/1000.;
        
        if(s.center.x() > maxx)
            maxx = s.center.x();
        if(s.center.y() > maxy)
            maxy = s.center.y();
        if(s.center.z() > maxz)
            maxz = s.center.z();
        if(s.center.x() < minx)
            minx = s.center.x();
        if(s.center.y() < miny)
            miny = s.center.y();
        if(s.center.z() < minz)
            minz = s.center.z();
        if(s.r > maxr)
            maxr = s.r;
        if(s.r < minr)
            minr = s.r;
        
        if (spheres_centers)
        {
            spheres_centers.write(reinterpret_cast<const char*>( &x ), sizeof( float ));
            spheres_centers.write(reinterpret_cast<const char*>( &z ), sizeof( float ));
            spheres_centers.write(reinterpret_cast<const char*>( &y ), sizeof( float ));
            spheres_centers.write(reinterpret_cast<const char*>( &r ), sizeof( float ));
        }
        //std::string center_str = std::to_string(s.center.x()/1000.) + " " + std::to_string(s.center.z()/1000.) + " " + std::to_string(s.center.y()/1000.) + " " + std::to_string(s.r/1000.);
        //const char* center_char = center_str.c_str();
        //if (spheres_centers) spheres_centers.write(center_char, sizeof(center_char));
        //const char* rho = std::to_string(s.rho).c_str();
        float rho = s.rho;
        if (densities) densities.write(reinterpret_cast<const char*>( &rho ), sizeof( float ));

        vx = s.speed.x();
        vy = s.speed.y();
        vz = s.speed.z();
        if (velocity)
        {
            velocity.write(reinterpret_cast<const char*>( &vx ), sizeof( float ));
            velocity.write(reinterpret_cast<const char*>( &vz ), sizeof( float ));
            velocity.write(reinterpret_cast<const char*>( &vy ), sizeof( float ));
        }
        
        chg = s.charge;
        if(charge)
            charge.write(reinterpret_cast<const char*>( &chg ), sizeof( float ));
        if(chg < minc)
            minc = chg;
        if(chg > maxc)
            maxc = chg;
        meanc += chg;
        countc++;
        
        int val = 2;
        if(sphtype)
            sphtype.write(reinterpret_cast<const char*>( &val ), sizeof( int ));
    }
    /*
    std::cerr << "FALLING SPHERES" << std::endl;
    std::cerr << "maxr = " << maxr << " minr = " << minr << std::endl;
    std::cerr << "maxx = " << maxx << " minx = " << minx << std::endl;
    std::cerr << "maxy = " << maxy << " miny = " << miny << std::endl;
    std::cerr << "maxz = " << maxz << " minz = " << minz << std::endl;
     */

    // falling spheres buffers
    minr = 100000.0f; maxr = 0.0f;
    minx = 100000.0f; maxx = -100000.0f; miny = 100000.0f; maxy = -100000.0f; minz = 100000.0f; maxz = -100000.0f;
    for (unsigned int k = 0; k<falling_spheres_buffers.size(); k++)
    {
        for (unsigned int i = 0; i<falling_spheres_buffers[k].size(); i++)
        {
            free_sphere_params s = falling_spheres_buffers[k][i];
            x = s.center.x()/1000.;
            y = s.center.y()/1000.;
            z = s.center.z()/1000.;
            r = s.r/1000.;
            
            if(s.center.x() > maxx)
                maxx = s.center.x();
            if(s.center.y() > maxy)
                maxy = s.center.y();
            if(s.center.z() > maxz)
                maxz = s.center.z();
            if(s.center.x() < minx)
                minx = s.center.x();
            if(s.center.y() < miny)
                miny = s.center.y();
            if(s.center.z() < minz)
                minz = s.center.z();
            if(s.r > maxr)
                maxr = s.r;
            if(s.r < minr)
                minr = s.r;
            
            if (spheres_centers)
            {
                spheres_centers.write(reinterpret_cast<const char*>( &x ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &z ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &y ), sizeof( float ));
                spheres_centers.write(reinterpret_cast<const char*>( &r ), sizeof( float ));
            }
            //std::string center_str = std::to_string(s.center.x()/1000.) + " " + std::to_string(s.center.z()/1000.) + " " + std::to_string(s.center.y()/1000.) + " " + std::to_string(s.r/1000.);
            //const char* center_char = center_str.c_str();
            //if (spheres_centers) spheres_centers.write(center_char, sizeof(center_char));
            //const char* rho = std::to_string(s.rho).c_str();
            float rho = s.rho;
            if (densities) densities.write(reinterpret_cast<const char*>( &rho ), sizeof( float ));
            
            vx = s.speed.x();
            vy = s.speed.y();
            vz = s.speed.z();
            if (velocity)
            {
                velocity.write(reinterpret_cast<const char*>( &vx ), sizeof( float ));
                velocity.write(reinterpret_cast<const char*>( &vz ), sizeof( float ));
                velocity.write(reinterpret_cast<const char*>( &vy ), sizeof( float ));
            }
            
            chg = 0.0f;
            if(charge)
                charge.write(reinterpret_cast<const char*>( &chg ), sizeof( float ));
            
            int val = 3;
            if(sphtype)
                sphtype.write(reinterpret_cast<const char*>( &val ), sizeof( int ));
        }
    }
    /*
    std::cerr << "SPHERE BUFFER" << std::endl;
    std::cerr << "maxr = " << maxr << " minr = " << minr << std::endl;
    std::cerr << "maxx = " << maxx << " minx = " << minx << std::endl;
    std::cerr << "maxy = " << maxy << " miny = " << miny << std::endl;
    std::cerr << "maxz = " << maxz << " minz = " << minz << std::endl;
     */
    
    // std::cerr << "CHARGE STATISTICS:" << std::endl;
    // std::cerr << "min = " << minc << " max = " << maxc << " mean = " << meanc / (float) countc << std::endl;

    spheres_centers.close();
    densities.close();
    velocity.close();
    charge.close();
    sphtype.close();
}


//------------------------------------------------------------
//------------------------- SETUP ----------------------------
//------------------------------------------------------------ */

void scene_model::setup_data(terrain_structure& terrain, std::shared_ptr<PWM::Model::planet> pl, unsigned int nb_layers)
{
    // begin with timer stopped
    //timer.stop();
    //replay = false;
    //frame_replay = 0;
    new_layer_delay = 0;
    total_layers_ejected = 0;
    nb_of_iterations = 0;
    frame_count = 0;
    run_time = 0;
    t_sum = 0;
    /*
    export_data = false;
    gui_param.display_smoke_layers = true;
    gui_param.display_free_spheres = false;
    gui_param.display_subspheres = false;
    gui_param.display_spheres_with_subspheres = false;
    gui_param.display_billboards = false;*/
    debug_mode = true;

    // seed
    // float seed = time(0);
    // srand(seed);
    /*if (seed_ofstream)
    {
        seed_ofstream.precision(10);
        seed_ofstream << seed;
    }*/

    //gui.show_frame_camera = false; std::cout << "replay becomes false 0" << std::endl;

    // camera setup
    /*scene.camera.apply_rotation(0, 0, 0, 1.2f);
    scene.camera.apply_scaling(100.0);
    scene.camera.apply_translation_in_screen_plane(0, -0.5);*/

    // Meshes setup
    /*
    layer_mesh = mesh_drawable( mesh_primitive_cylinder(0.1f, {0,0,0}, {0,0,0.01}));
    layer_mesh.shader = shaders["mesh"];
    layer_mesh.uniform.color = {0,0.5,1};

    mesh cyl = vcl::mesh_primitive_cylinder(2.5f, {0,0,1.5}, {0,0,-1.5}, 30, 30);
    mesh d1 = vcl::mesh_primitive_disc(2.5f, {0,0,1.5});
    mesh d2 = vcl::mesh_primitive_disc(2.5f, {0,0,-1.5});
    mesh t = cyl;t.push_back(d1); t.push_back(d2);

    generic_sphere_mesh = vcl::mesh_primitive_sphere();
    generic_sphere_mesh.texture_id = scene.texture_white;
    //generic_torus_mesh = vcl::mesh_primitive_torus(1.,1.,{0,0,0}, {0,0,-1});
    generic_torus_mesh = t;
    generic_torus_mesh.uniform.color = {1,0.5,0};
    generic_torus_mesh.shader = shaders["mesh"];
    generic_torus_mesh.texture_id = scene.texture_white;
    generic_torus_mesh.uniform.color_alpha = 0.6f;
    texture_smoke_id = create_texture_gpu( image_load_png("../scenes/sources/smoke/images/texture_panache.png") );

    sphere = mesh_drawable( mesh_primitive_sphere(0.1f));
    sphere.shader = shaders["mesh"];
    sphere.uniform.color = {0,0.5,1};
    sphere.texture_id = scene.texture_white;

    subspheres = vcl::mesh_primitive_sphere(1.0, {0,0,0}, 10 ,20);
    subspheres.texture_id = scene.texture_white;
    subspheres.uniform.color = {0.6,0.5,0.5};
    subspheres.uniform.shading.diffuse = 0.8f;
    subspheres.uniform.shading.specular = 0.0f;

    GLuint const texture_billboard = create_texture_gpu(image_load_png("../scenes/sources/smoke/terrains/smoke2.png"));
    quad = mesh_drawable(mesh_primitive_quad({-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0}));
    quad.texture_id = texture_billboard;
    quad.uniform.shading.ambiant = 1.0;
    quad.uniform.shading.diffuse = 0.0;
    quad.uniform.shading.specular = 0.0;

    auto circle = vcl::curve_primitve_circle(30, 1.0, {0,0,0}, {0,0,1});
    sphere_circle = curve_drawable(circle);
    sphere_circle.shader = shaders["curve"];
    sphere_circle.uniform.color = {1,0,0};

    //sampling subpheres
    {
        int N = 60;
        for (int k = 0; k < N; ++k)
        {
            //uniform sampling on sphere
            float theta = 2*3.14f*vcl::rand_interval();
            float phi   = std::acos(1-2.0f*vcl::rand_interval());


            float x = std::sin(phi)*std::cos(theta);
            float y = std::sin(phi)*std::sin(theta);
            float z = std::cos(phi);

            vec3 p = {x,y,z};
            bool add = true;
            for (int k2 = 0; add==true && k2 < k; ++k2)
                if(norm(p-samples_subspheres[k2])<0.18f)
                    add=false;
            samples_subspheres.push_back({x,y,z});
        }
    }

    {
        mesh m0 = vcl::mesh_primitive_sphere(1.0, {0,0,0}, 5 ,5);
        mesh m1 = vcl::mesh_primitive_sphere(1.0, {0,0,0}, 8 ,8);
        mesh m2 = vcl::mesh_primitive_sphere(1.0, {0,0,0}, 10 , 10);

        int N = 60;
        for (int k = 0; k < N; ++k)
        {
            //uniform sampling on sphere
            float theta = 2*3.14f*vcl::rand_interval();
            float phi   = std::acos(1-2.0f*vcl::rand_interval());
            float r = vcl::rand_interval(0.8f,1.0f);

            float x = r*std::sin(phi)*std::cos(theta);
            float y = r*std::sin(phi)*std::sin(theta);
            float z = r*std::cos(phi);

            vec3 p = {x,y,z};
            bool add = true;
            for (int k2 = 0; add==true && k2 < k; ++k2)
                if(norm(p-samples_subspheres[k2])<0.18f)
                    add=false;
            samples_subspheres.push_back({x,y,z});
        }

        mesh m;
        m.push_back(m0);
        for (int sub = 0; sub < samples_subspheres.size(); ++sub) {
            mesh temp = m1;
            float r = vcl::rand_interval(0.18f,0.2f);

            // subspheres
            for (int k = 0; k < temp.position.size(); ++k)
            {
                vec3 p = r*temp.position[k] + samples_subspheres[sub];

                vec3 n0 = temp.normal[k];
                vec3 n1 = normalize(p);

                float d = norm(p);
                float alpha = 0.0;
                if(d>1.0f && d<1.2f)
                    alpha = (d-1.0f)/0.2f;
                if(d>1.2f)
                    alpha = 1.0f;

                vec3 n = (1-alpha)*n1 + alpha*n0; // hack normals
                temp.normal[k] = n;
                temp.position[k] = p;
            }

            m.push_back(temp);
        }

        subspheres_display = mesh_drawable(m) ;

        subspheres_display.texture_id = scene.texture_white;
        subspheres_display.uniform.color = {0.6,0.6,0.55};
        subspheres_display.uniform.shading.ambiant = 0.7f;
        subspheres_display.uniform.shading.diffuse = 0.3f;
        subspheres_display.uniform.shading.specular = 0.0f;
    }
    */

    // Terrain setup
    terrain_struct = terrain;
    planetos = pl;
    /* mesh mesh_terrain = mesh_load_file_obj("../scenes/sources/smoke/terrains/sthelens_detailed_sub5.obj");
    // Fill std::vector<Eigen::Vector3f> for positions and normals)
    std::vector<Eigen::Vector3f> positions;
    std::vector<Eigen::Vector3f> normals;
    for (unsigned int i = 0; i < mesh_terrain.position.size(); i++)
    {
        Eigen::Vector3f pos = mesh_terrain.position[i];
        positions.push_back(pos);

        Eigen::Vector3f nor = mesh_terrain.normal[i];
        normals.push_back(nor);
    }
    //terrain = mesh_drawable(mesh_terrain);
    //terrain.shader = shaders["mesh"];
    //terrain.uniform.color = {0.36f, 0.32f, 0.24f};
    //terrain.uniform.shading.specular = 0.0f;
    Eigen::Matrix3f rotation = (Eigen::AngleAxisf(3.14f/2.0f, Eigen::Vector3f(1.0f, 0, 0))).toRotationMatrix();
    float scaling = 100.f;
    Eigen::Vector3f translation(37.f,68.f,-25.f);
    fill_height_field(positions, normals, rotation, scaling, translation); */

    //terrain_display = terrain;
    //terrain_display.uniform.transform.scaling = 1./1.;
    //terrain_display.texture_id = create_texture_gpu(image_load_png("../scenes/sources/smoke/terrains/terrain.png"));
    //terrain_display.uniform.color = {1,1,1};

    // Params setup
    is_wind = false;
    linear_wind_base = 15.;
    for(unsigned int i = 0; i< nb_layers; i++)
    {
        wind_altitudes.push_back(i*4000);
        winds.push_back(wind_structure(0,0));
        temperatures.push_back(compute_atm_temperature(i*4000.));
        moistures.push_back(0);
        pressures.push_back(compute_atm_pressure(i*4000.));
    }

    stagnation_speed = 50;

    // coeff init
    air_incorporation_coeff = 5.;

    subspheres_number = 300;
    subsubspheres_number = 0;

    //reset terrain heights to match PWM with min = 0 m;
    float hOffset = std::numeric_limits<float>::max();
    for (int i = 0; i < terrain_struct.height_field.size(); ++i){
        for (int j = 0; j < terrain_struct.height_field[i].size(); ++j){
            if (terrain_struct.height_field[i][j] < hOffset)
                hOffset = terrain_struct.height_field[i][j];
        }
    }
    for (int i = 0; i < terrain_struct.height_field.size(); ++i){
        for (int j = 0; j < terrain_struct.height_field[i].size(); ++j){
            terrain_struct.height_field[i][j] = terrain_struct.height_field[i][j] - hOffset;
        }
    }

    // Parameters : to be chosen by user
    T_0 = 1273.; // initial temp (K)
    theta_0 = 0.; // initial angle (rad)
    U_0 = 200.; // initial speed (m.s-1)
    n_0 = 0.03; // initial gas mass fraction
    z_0 = 0 - hOffset; // initial altitude (m) //PCP Fix to read altitude from known height.
    r_0 = 150.; // initial radius (m) // JG change from 100 to 150
    c_0 = 3.0; // initial charge
    rho_0 = 200.;
    
    iceChargeMax = 50.0f; // 50 g/m^3 of water as the maximum threshold
    ashChargeColdMax = 5.0f; // ash density threshold above -20C isotherm
    ashChargeHotMax = 100.0f; // ash density threshold below -20C isotherm
    convectChargeMax = 20.0f; // 5 m/s convection speed threshold
    coldChargeCoeff = 0.5f; // proportion of contribution to charge is 10X for cold zone
    hotChargeCoeff = 0.1f;
    velChargeMax = 400.0f; // velocity threshold in jet phase
    moistureChargeMax = 51.0f; // water vapour amount for 100% relative humidity at 40C

    // Parameters : constants
    g = planetos->getGravitationalAcceleration(); // (m.s-2)
}

/*
void scene_model::display(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    draw(terrain_display, scene.camera, shaders["mesh"]);
    //draw(terrain, scene.camera, shaders["wireframe"]);

    float ratio = 100.0;

    // Display torus
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for (unsigned int i=0; i<smoke_layers.size(); i++)
    {
        smoke_layer lay = smoke_layers[i];

        generic_torus_mesh.uniform.transform.scaling = lay.r/ratio;
        generic_torus_mesh.uniform.transform.translation = vec3(lay.center.x/ratio, lay.center.y/ratio, lay.center.z/ratio);
        generic_torus_mesh.uniform.transform.rotation = rotation_from_axis_angle_mat3(lay.theta_axis, lay.theta-3.14/2.0);
        if(gui_param.display_smoke_layers) draw(generic_torus_mesh, scene.camera);
    }

//    glBindTexture(GL_TEXTURE_2D, texture_smoke_id);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

    // billboards
    if(gui_param.display_billboards)
    {
        glDepthMask(false);
        for (unsigned int j = 0; j<free_spheres.size(); j++)
        {
            mat3 const R = rotation_from_axis_angle_mat3(free_spheres[j].rotation_axis, free_spheres[j].current_angle);
            float new_scaling = free_spheres[j].r/ratio;
            //if (j==0) std::cout << new_scaling << std::endl;
            vec3 new_translation = free_spheres[j].center/ratio;
            generic_sphere_mesh.uniform.transform.translation = new_translation;
            generic_sphere_mesh.uniform.transform.scaling = new_scaling;
            generic_sphere_mesh.uniform.transform.rotation = R;
            generic_sphere_mesh.uniform.color = {1,1,1};

            float var = vcl::perlin(j,2);

            quad.uniform.transform.rotation = rotation_from_axis_angle_mat3(scene.camera.orientation.col(2), free_spheres[j].current_angle * dot(free_spheres[j].rotation_axis,scene.camera.orientation.col(2)) * 1.5f *(1+0.3*var)   + 2.2145*j*j) * scene.camera.orientation;
            quad.uniform.transform.translation = new_translation;
            quad.uniform.transform.scaling = new_scaling*1.3;
            quad.uniform.color_alpha = 0.8+0.3f*(2*var-1.0f);
            draw(quad, scene.camera, shaders["mesh"]);
        }
        glDepthMask(true);
    }


    // free + stagnation spheres display
    if(gui_param.display_free_spheres)
    {
        for (unsigned int j = 0; j<free_spheres.size(); j++)
        {
            //if (!free_spheres[j].falling)
            if(true)
            {
                mat3 const R = rotation_from_axis_angle_mat3(free_spheres[j].rotation_axis, free_spheres[j].current_angle);
                float r = free_spheres[j].r/ratio;
                vec3 t = free_spheres[j].center/ratio;
                float rho = free_spheres[j].rho;
                float disp_rho = 1. - rho;
                if (disp_rho < 0) disp_rho = 0.;
                disp_rho = 1.;

                generic_sphere_mesh.uniform.transform.translation = t;
                generic_sphere_mesh.uniform.transform.scaling = r;
                generic_sphere_mesh.uniform.transform.rotation = R;
                generic_sphere_mesh.uniform.color = {disp_rho,disp_rho,disp_rho};
                if(gui_param.display_free_spheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
            }
        }
    }

    // spheres+subspheres display (lighter)
    if(gui_param.display_spheres_with_subspheres)
    {
        for (unsigned int j = 0; j<free_spheres.size(); j++)
        {
            //if (!free_spheres[j].falling)
            if (true)
            {
                mat3 const R = rotation_from_axis_angle_mat3(free_spheres[j].rotation_axis, free_spheres[j].current_angle);
                float r = free_spheres[j].r/ratio;
                vec3 t = free_spheres[j].center/ratio;
                float rho = free_spheres[j].rho;
                float disp_rho = 1. - rho;
                if (disp_rho < 0) disp_rho = 0.;
                disp_rho = 1.;

                subspheres_display.uniform.transform.translation = t;
                subspheres_display.uniform.transform.scaling = r;
                subspheres_display.uniform.transform.rotation = R;
                subspheres_display.uniform.color = {disp_rho,disp_rho,disp_rho};

                draw(subspheres_display, scene.camera, shaders["mesh"]);
            }
        }
    }

    // subspheres display
    if(gui_param.display_subspheres)
    {
        for (unsigned int j = 0; j<s2_spheres.size(); j++)
        {
            if (!free_spheres[s2_spheres[j].parent_id].falling)
            {
                float r = s2_spheres[j].r/ratio;
                vec3 t = s2_spheres[j].center/ratio;
                float rho = free_spheres[s2_spheres[j].parent_id].rho;
                float disp_rho = 1. - rho;
                if (disp_rho < 0) disp_rho = 0.;
                disp_rho = 1.;

                generic_sphere_mesh.uniform.transform.translation = t;
                generic_sphere_mesh.uniform.transform.scaling = r;
                generic_sphere_mesh.uniform.color = {disp_rho,disp_rho,disp_rho};
                if(gui_param.display_subspheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
            }
        }
    }

    // falling spheres display
    for (unsigned int j = 0; j<falling_spheres.size(); j++)
    {
        float new_scaling = falling_spheres[j].r/ratio;
        vec3 new_translation = {falling_spheres[j].center.x/ratio, falling_spheres[j].center.y/ratio, falling_spheres[j].center.z/ratio};
        generic_sphere_mesh.uniform.transform.translation = new_translation;
        generic_sphere_mesh.uniform.transform.scaling = new_scaling;
        generic_sphere_mesh.uniform.transform.rotation = mat3::identity();
        if (falling_spheres[j].falling_under_atm_rho) generic_sphere_mesh.uniform.color = {1,0,0};
        else generic_sphere_mesh.uniform.color = {1,1,1};
        if(gui_param.display_free_spheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
    }
    // buffer falling spheres display
    for (unsigned int k = 0; k<falling_spheres_buffers.size(); k++)
    {
        for (unsigned int j = 0; j<falling_spheres_buffers[k].size(); j++)
        {
            float new_scaling = falling_spheres_buffers[k][j].r/ratio;
            vec3 new_translation = {falling_spheres_buffers[k][j].center.x/ratio, falling_spheres_buffers[k][j].center.y/ratio, falling_spheres_buffers[k][j].center.z/ratio};
            generic_sphere_mesh.uniform.transform.translation = new_translation;
            generic_sphere_mesh.uniform.transform.scaling = new_scaling;
            generic_sphere_mesh.uniform.transform.rotation = mat3::identity();
            if (falling_spheres_buffers[k][j].falling_under_atm_rho) generic_sphere_mesh.uniform.color = {1,0,0};
            else generic_sphere_mesh.uniform.color = {1,1,1};
            if(gui_param.display_free_spheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
        }
    }

    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

}

void scene_model::display_replay(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui, size_t frame)
{
    draw(terrain_display, scene.camera, shaders["mesh"]);

    float ratio = 100.0;

    // Display torus
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for (unsigned int i=0; i<smoke_layers_frames[frame].size(); i++)
    {
        smoke_layer lay = smoke_layers_frames[frame][i];

        generic_torus_mesh.uniform.transform.scaling = lay.r/ratio;
        generic_torus_mesh.uniform.transform.translation = vec3(lay.center.x/ratio, lay.center.y/ratio, lay.center.z/ratio);
        generic_torus_mesh.uniform.transform.rotation = rotation_from_axis_angle_mat3(lay.theta_axis, lay.theta-3.14/2.0);
        if(gui_param.display_smoke_layers) draw(generic_torus_mesh, scene.camera);
    }


//    glBindTexture(GL_TEXTURE_2D, texture_smoke_id);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

    // billboards
    if(gui_param.display_billboards)
    {
        glDepthMask(false);
        for (unsigned int j = 0; j<free_spheres_frames[frame].size(); j++)
        {
            mat3 const R = rotation_from_axis_angle_mat3(free_spheres_frames[frame][j].rotation_axis, free_spheres_frames[frame][j].current_angle);
            float new_scaling = free_spheres_frames[frame][j].r/ratio;
            vec3 new_translation = free_spheres_frames[frame][j].center/ratio;
            generic_sphere_mesh.uniform.transform.translation = new_translation;
            generic_sphere_mesh.uniform.transform.scaling = new_scaling;
            generic_sphere_mesh.uniform.transform.rotation = R;
            generic_sphere_mesh.uniform.color = {1,1,1};

            float var = vcl::perlin(j,2);

            quad.uniform.transform.rotation = rotation_from_axis_angle_mat3(scene.camera.orientation.col(2), free_spheres_frames[frame][j].current_angle * dot(free_spheres_frames[frame][j].rotation_axis,scene.camera.orientation.col(2)) * 1.5f *(1+0.3*var)   + 2.2145*j*j) * scene.camera.orientation;
            quad.uniform.transform.translation = new_translation;
            quad.uniform.transform.scaling = new_scaling;
            quad.uniform.color_alpha = 0.8+0.3f*(2*var-1.0f);
            draw(quad, scene.camera, shaders["mesh"]);
        }
        glDepthMask(true);
    }

    // free + stagnation spheres display NEW
    for (unsigned int j = 0; j<free_spheres_frames[frame].size(); j++)
    {
        if (!free_spheres_frames[frame][j].falling)
        {
            mat3 const R = rotation_from_axis_angle_mat3(free_spheres_frames[frame][j].rotation_axis, free_spheres_frames[frame][j].current_angle);
            float r = free_spheres_frames[frame][j].r/ratio;
            vec3 t = free_spheres_frames[frame][j].center/ratio;
            float rho = free_spheres_frames[frame][j].rho;
            float disp_rho = 1. - rho;
            if (disp_rho < 0) disp_rho = 0.;

            generic_sphere_mesh.uniform.transform.translation = t;
            generic_sphere_mesh.uniform.transform.scaling = r;
            generic_sphere_mesh.uniform.transform.rotation = R;
            generic_sphere_mesh.uniform.color = {disp_rho,disp_rho,disp_rho};
            if(gui_param.display_free_spheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
        }
    }

    // spheres+subspheres display (lighter)
    if(gui_param.display_spheres_with_subspheres)
    {
        for (unsigned int j = 0; j<free_spheres_frames[frame].size(); j++)
        {
            if (!free_spheres_frames[frame][j].falling)
            {
                mat3 const R = rotation_from_axis_angle_mat3(free_spheres_frames[frame][j].rotation_axis, free_spheres_frames[frame][j].current_angle);
                float r = free_spheres_frames[frame][j].r/ratio;
                vec3 t = free_spheres_frames[frame][j].center/ratio;
                float rho = free_spheres_frames[frame][j].rho;
                float disp_rho = 1. - rho;
                if (disp_rho < 0) disp_rho = 0.;

                subspheres_display.uniform.transform.translation = t;
                subspheres_display.uniform.transform.scaling = r;
                subspheres_display.uniform.transform.rotation = R;
                subspheres_display.uniform.color = {disp_rho,disp_rho,disp_rho};

                draw(subspheres_display, scene.camera, shaders["mesh"]);
            }
        }
    }

    // subspheres display NEW
    if(gui_param.display_subspheres)
    {
        for (unsigned int j = 0; j<s2_spheres_frames[frame].size(); j++)
        {
            if (!free_spheres[s2_spheres_frames[frame][j].parent_id].falling)
            {
                float r = free_spheres_frames[frame][j].r/ratio;
                vec3 t = free_spheres_frames[frame][j].center/ratio;
                float rho = free_spheres_frames[frame][j].rho;
                float disp_rho = 1. - rho;
                if (disp_rho < 0) disp_rho = 0.;

                generic_sphere_mesh.uniform.transform.translation = t;
                generic_sphere_mesh.uniform.transform.scaling = r;
                generic_sphere_mesh.uniform.color = {disp_rho,disp_rho,disp_rho};
                if(gui_param.display_subspheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
            }
        }
    }

    // falling spheres display
    for (unsigned int j = 0; j<falling_spheres_frames[frame].size(); j++)
    {
        float new_scaling = falling_spheres_frames[frame][j].r/ratio;
        vec3 new_translation = {falling_spheres_frames[frame][j].center.x/ratio, falling_spheres_frames[frame][j].center.y/ratio, falling_spheres_frames[frame][j].center.z/ratio};
        generic_sphere_mesh.uniform.transform.translation = new_translation;
        generic_sphere_mesh.uniform.transform.scaling = new_scaling;
        generic_sphere_mesh.uniform.transform.rotation = mat3::identity();
        if (falling_spheres_frames[frame][j].falling_under_atm_rho) generic_sphere_mesh.uniform.color = {1,0,0};
        else generic_sphere_mesh.uniform.color = {1,1,1};
        if(gui_param.display_free_spheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
    }
    // buffer falling spheres display
    for (unsigned int k = 0; k<falling_spheres_buffers_frames[frame].size(); k++)
    {
        for (unsigned int j = 0; j<falling_spheres_buffers_frames[frame][k].size(); j++)
        {
            float new_scaling = falling_spheres_buffers_frames[frame][k][j].r/ratio;
            vec3 new_translation = {falling_spheres_buffers_frames[frame][k][j].center.x/ratio, falling_spheres_buffers_frames[frame][k][j].center.y/ratio, falling_spheres_buffers_frames[frame][k][j].center.z/ratio};
            generic_sphere_mesh.uniform.transform.translation = new_translation;
            generic_sphere_mesh.uniform.transform.scaling = new_scaling;
            generic_sphere_mesh.uniform.transform.rotation = mat3::identity();
            if (falling_spheres_buffers_frames[frame][k][j].falling_under_atm_rho) generic_sphere_mesh.uniform.color = {1,0,0};
            else generic_sphere_mesh.uniform.color = {1,1,1};
            if(gui_param.display_free_spheres) draw(generic_sphere_mesh, scene.camera, shaders["mesh"]);
        }
    }

    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
}


void scene_model::set_gui()
{
    // Can set the speed of the animation
    float scale_min = 0.05f;
    float scale_max = 2.0f;
    ImGui::SliderScalar("Time scale", ImGuiDataType_Float, &timer.scale, &scale_min, &scale_max, "%.2f s");

    // Parameters

    ImGui::Checkbox("Display torus layers", &gui_param.display_smoke_layers);
    ImGui::Checkbox("Display free spheres", &gui_param.display_free_spheres);
    //ImGui::Checkbox("Display subspheres", &gui_param.display_subspheres);
    ImGui::Checkbox("Display spheres with subspheres", &gui_param.display_spheres_with_subspheres);
    ImGui::Checkbox("Display billboards", &gui_param.display_billboards);

    unsigned int spheres_min=0, spheres_max=500;
    ImGui::SliderScalar("Number of subspheres", ImGuiDataType_S32, &subspheres_number, &spheres_min, &spheres_max);
    ImGui::SliderScalar("Number of subsubspheres", ImGuiDataType_S32, &subsubspheres_number, &spheres_min, &spheres_max);

    // Coeffs

    //float air_inc_min = 0.5, air_inc_max = 10.;
    //ImGui::SliderScalar("Air incorporation coefficient", ImGuiDataType_Float, &air_incorporation_coeff, &air_inc_min, &air_inc_max, "%.2f");

    // Initial conditions

    float initial_speed_min = 0., initial_speed_max = 200.;
    ImGui::SliderScalar("Initial plume speed", ImGuiDataType_Float, &U_0, &initial_speed_min, &initial_speed_max, "%.2f m/s");
    float initial_density_min = 150., initial_density_max = 250.;
    ImGui::SliderScalar("Initial plume density", ImGuiDataType_Float, &rho_0, &initial_density_min, &initial_density_max, "%.2f kg/m3");
    float vent_ray_min = 50., vent_ray_max = 200.;
    ImGui::SliderScalar("Vent radius", ImGuiDataType_Float, &r_0, &vent_ray_min, &vent_ray_max, "%.2f m");
    float vent_altitude_min = 0., vent_altitude_max = 8000.;
    ImGui::SliderScalar("Vent altitude", ImGuiDataType_Float, &z_0, &vent_altitude_min, &vent_altitude_max, "%.2f m");


    // Wind presets
    if (ImGui::Button("No wind"))
    {
        is_wind = false;
        for(unsigned int i = 0; i<winds.size(); i++)
        {
            winds[i].intensity = 0;
            winds[i].angle = 0;
            winds[i].wind_vector = vec3(1,0,0);
        }
    }
    if (ImGui::Button("Linear Wind"))
    {
        is_wind = true;
        for(unsigned int i = 0; i<winds.size(); i++)
        {
            winds[i].intensity = i*linear_wind_base;
            if (i>3) winds[i].intensity = 3*linear_wind_base;
            if (winds[i].intensity == 0) winds[i].intensity = 1;
            winds[i].angle = 0;
            winds[i].wind_vector = winds[i].intensity * vec3(cos(winds[i].angle), sin(winds[i].angle),0);
        }
    }

    // Wind

    float lin_windbase_min = 0., lin_windbase_max = 35.;
    if (ImGui::SliderScalar("Linear wind base speed", ImGuiDataType_Float, &linear_wind_base, &lin_windbase_min, &lin_windbase_max, "%1.f m/s"))
    {
        if (is_wind)
        {
            for(unsigned int i = 0; i<winds.size(); i++)
            {
                winds[i].intensity = i*linear_wind_base;
                if (i>3) winds[i].intensity = 3*linear_wind_base;
                if (winds[i].intensity == 0) winds[i].intensity = 1;
                winds[i].angle = 0;
                winds[i].wind_vector = winds[i].intensity * vec3(cos(winds[i].angle), sin(winds[i].angle),0);
            }
        }
    }

    int wind_min = 0;
    int wind_max = 300;
    int angle_min = 0, angle_max = 360;
    for (int i = wind_altitudes.size()-1; i>=0; i--)
    {
        std::string alti = "Intensity (" + std::to_string(wind_altitudes[i]) + "m)";
        std::string angl = "Angle (" + std::to_string(wind_altitudes[i]) + "m)";
        if (ImGui::SliderScalar(alti.c_str(), ImGuiDataType_S32, &winds[i].intensity, &wind_min, &wind_max))
            winds[i].wind_vector = winds[i].intensity * vec3(cos(winds[i].angle), sin(winds[i].angle),0);
        if (ImGui::SliderScalar(angl.c_str(), ImGuiDataType_S32, &winds[i].angle, &angle_min, &angle_max))
            winds[i].wind_vector = winds[i].intensity * vec3(cos(winds[i].angle), sin(winds[i].angle),0);
    }

    // Start and stop animation
    if (ImGui::Button("Stop"))
        timer.stop();
    if (ImGui::Button("Start"))
        timer.start();
    if (ImGui::Button("Stop and reset"))
    {
        timer.stop();
        smoke_layers.clear();
        free_spheres.clear();
        falling_spheres.clear();
        stagnate_spheres.clear();
        falling_spheres_buffers.clear();
        frame_count = 0;

        smoke_layers_frames.clear();
        free_spheres_frames.clear();
        s2_spheres_frames.clear();
        stagnate_spheres_frames.clear();
        falling_spheres_frames.clear();
        falling_spheres_buffers_frames.clear();
        frame_replay = 0;
        replay = false;
        export_data = false;
    }
    if (ImGui::Button("Replay"))
    {
        if (!export_data)
        {
            timer.stop();
            frame_replay = 0;
            replay = true; std::cout << "replay becomes true" << std::endl;
        }
    }
    if (ImGui::Button("Stop Replay"))
    {
        timer.start();
        replay = false; std::cout << "replay becomes false" << std::endl;
    }
    if (ImGui::Button("Export data during simulation (slow, no replay, activate before starting simulation)"))
    {
        export_data = true;
        replay = false;
    }

}

*/

