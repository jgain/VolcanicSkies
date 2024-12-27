#ifndef CF_SMOKE_HPP
#define CF_SMOKE_HPP

#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include "flatStaggeredGrid.h"
#include <iostream>
#include <fstream>
#include "planet.h"
#include "smokeLayer.hpp"
#include <string>
#include "terrain_structure.hpp"
#include "world.h"
#include <thread>

typedef double valType;
typedef std::string valType2;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType2> dsSType;

struct wind_structure
{
    float intensity;
    float angle;
    Eigen::Vector3f wind_vector; // horizontal

    wind_structure() : intensity(0), angle(0), wind_vector(1,0,0) {}
    wind_structure(float intensity, float angle) : intensity(intensity), angle(angle)
    {
        wind_vector = intensity * Eigen::Vector3f(cos(angle), sin(angle), 0);
    }
    wind_structure(float x, float y, bool are_coord)
    {
        intensity = sqrt(x*x+y*y);
        if (x < 0.001 && x > -0.001) angle = 3.14/2.;
        else angle = atan(y/x);
        wind_vector = Eigen::Vector3f(x,y,0);
    }
};

// User parameters available in the GUI
/*
struct gui_parameters
{
    bool display_smoke_layers;
    bool display_free_spheres;
    bool display_subspheres;
    bool display_spheres_with_subspheres;
    bool display_billboards;
};*/


struct scene_model
{
    unsigned int frame_count;
    float run_time;
    bool debug_mode;
    float dt;
    float t_sum;
    bool replay;
    size_t frame_replay;
    bool export_data;
    
    // for tracking charge distribution
    float chrghotsum,  chrgcoldsum;
    int chrghotcnt, chrgcoldcnt;

    std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>> atmoModel;

    // Trackers
    float new_layer_delay;
    unsigned int total_layers_ejected;
    unsigned int nb_of_iterations;
    unsigned int last_ppe_layer_idx;

    // Meshes
    /*
    vcl::mesh_drawable generic_torus_mesh;
    vcl::mesh_drawable generic_sphere_mesh;
    vcl::mesh_drawable layer_mesh;
    vcl::mesh_drawable terrain;
    vcl::mesh_drawable terrain_display;
    vcl::mesh_drawable sphere;
    vcl::mesh_drawable quad;
    vcl::curve_drawable sphere_circle;
    GLuint texture_smoke_id;

    std::vector<vcl::vec3> samples_subspheres;
    vcl::mesh_drawable subspheres;
    vcl::mesh_drawable subspheres_display;*/

    // Parameters : to be chosen by user
    float T_0; // initial temp
    float theta_0; // initial angle
    float U_0; // initial speed
    float n_0; // initial gas mass fraction
    float z_0; // initial altitude
    float r_0; // initial radius
    float rho_0; // initial density
    float c_0; // initial charge
    float air_incorporation_coeff;
    float stagnation_speed;
    Eigen::Vector3f plume_loc;
    
    // charge formula coefficients
    float iceChargeMax, ashChargeColdMax, ashChargeHotMax, convectChargeMax, coldChargeCoeff, hotChargeCoeff, velChargeMax, moistureChargeMax;

    unsigned int subspheres_number;
    unsigned int subsubspheres_number;

    std::vector<int> wind_altitudes;
    std::vector<wind_structure> winds;
    std::vector<float> temperatures;
    std::vector<float> pressures;
    std::vector<float> moistures;
    float linear_wind_base;
    bool is_wind;

    // Parameters : constants
    float g;

    // Data structures
    std::vector<smoke_layer> smoke_layers;
    std::vector<free_sphere_params> free_spheres;
    std::vector<subsphere_params> s2_spheres;
    std::vector<subsphere_params> s3_spheres;
    std::vector<free_sphere_params> stagnate_spheres;
    std::vector<free_sphere_params> falling_spheres;
    std::vector< std::vector<free_sphere_params> > falling_spheres_buffers;
    terrain_structure terrain_struct;
    std::shared_ptr<PWM::Model::planet> planetos;

    // For replay feature
    /*
    std::vector< std::vector<smoke_layer> > smoke_layers_frames;
    std::vector< std::vector<free_sphere_params> > free_spheres_frames;
    std::vector< std::vector<subsphere_params> > s2_spheres_frames;
    std::vector< std::vector<free_sphere_params> > stagnate_spheres_frames;
    std::vector< std::vector<free_sphere_params> > falling_spheres_frames;
    std::vector< std::vector< std::vector<free_sphere_params> > > falling_spheres_buffers_frames;*/

    // Export structures
    std::string const altitude_file = "../output/altitude.txt";
    std::string const plumex_file = "../output/plumex_file.txt";
    std::string const speed_file = "../output/speed.txt";
    std::string const ray_file = "../output/ray.txt";
    std::string const smoke_rho_file = "../output/smoke_rho.txt";
    std::string const atm_rho_file = "../output/atm_rho.txt";
    std::string const temp_file = "../output/temp.txt";
    std::string const moist_file = "../output/moist.txt";
    std::string const time_file = "../output/time.txt";
    std::ofstream altitude = std::ofstream(altitude_file.c_str());
    std::ofstream plumex = std::ofstream(plumex_file.c_str());
    std::ofstream speed = std::ofstream(speed_file.c_str());
    std::ofstream ray = std::ofstream(ray_file.c_str());
    std::ofstream smoke_rho = std::ofstream(smoke_rho_file.c_str());
    std::ofstream atm_rho = std::ofstream(atm_rho_file.c_str());
    std::ofstream temp = std::ofstream(temp_file.c_str());
    std::ofstream moist = std::ofstream(moist_file.c_str());
    std::ofstream time = std::ofstream(time_file.c_str());
    // Seed export
    std::string const seed_file = "../output/seed.txt";
    std::ofstream seed_ofstream = std::ofstream(seed_file.c_str());

    // Test T and rho
    std::string const theoric_temp_file = "../output/theoric_temp.txt";
    std::string const actual_temp_file = "../output/actual_temp.txt";
    std::string const theoric_rho_file = "../output/theoric_rho.txt";
    std::string const actual_rho_file = "../output/actual_rho.txt";
    std::string const theoric_pres_file = "../output/theoric_pres.txt";
    std::string const actual_pres_file = "../output/actual_pres.txt";
    std::ofstream theoric_temp = std::ofstream(theoric_temp_file.c_str());
    std::ofstream actual_temp = std::ofstream(actual_temp_file.c_str());
    std::ofstream theoric_rho = std::ofstream(theoric_rho_file.c_str());
    std::ofstream actual_rho = std::ofstream(actual_rho_file.c_str());
    std::ofstream theoric_pres = std::ofstream(theoric_pres_file.c_str());
    std::ofstream actual_pres = std::ofstream(actual_pres_file.c_str());


    // General functions
    void setup_data(terrain_structure& terrainData, std::shared_ptr<PWM::Model::planet> pl, unsigned int nb_layers);
    void frame_draw(float deltaT);
    //void display(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    //void display_replay(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui, size_t frame);

    // Smoke layer computation
    float compute_atm_temperature(float height);
    float compute_atm_density(float height);
    float compute_atm_pressure(float height);
    Eigen::Vector3f compute_wind_vector(float height);
    float compute_temperature_at_h(float height);
    float compute_moisture_at_h(float height);
    float compute_density_at_h(float height);
    float compute_pressure_at_h(float height);
    void add_smoke_layer(float v, float d, float c, float r, Eigen::Vector3f position, bool secondary_plume);
    void edit_smoke_layer_properties(unsigned int i, float& d_mass);
    void apply_forces_to_smoke_layer(unsigned int i, float d_mass);
    void sedimentation(unsigned int i, float& d_mass);
    void pyroclastic_flow_computation_step(unsigned int i);
    void complete_plume_layer_properties_update(unsigned int i);
    void smoke_layer_update(unsigned int i);

    // Pyroclastic flow : falling spheres
    float field_height_at(float x, float y);
    Eigen::Vector3f field_normal_at(float x, float y);
    void sphere_ground_collision(free_sphere_params& sphere, int idx, unsigned int frame_nb);
    void ground_falling_sphere_update(free_sphere_params& sphere, int idx, unsigned int frame_nb);
    void secondary_columns_creation();
    void falling_spheres_update(unsigned int frame_nb);

    // Free spheres
    void add_free_sphere(unsigned int i, float angle, float size_fac);
    void add_free_spheres_for_one_layer(unsigned int i);
    void subdivide_and_make_falling(unsigned int i);
    void update_free_spheres();

    // Stagnation
    void update_stagnation_spheres();
    Eigen::Vector2f sampleVel(Eigen::Vector3f& loc);
    
    // Skirts
    void boxSpheres();
    void boundSpheres(float zone); // minz and maxz for a sphere trajectory in the last timestep
    int calcSubsteps(); // subdivision of main timestep to avoid sampling errors for skirts
    void trajectoryInterp(Eigen::Vector3f& loc, Eigen::Vector3f& vel, Eigen::Vector3f& subloc, Eigen::Vector3f& subvel, float t); // subsample sphere trajectory
    float calcFallOff(Eigen::Vector3f& loc, Eigen::Vector3f& c, float r, float hzone, float vzone);
    float calcSkirtUplift(Eigen::Vector3f& loc, float deltaT, float hzone, float vzone, int substeps);
    
    // Charge
    
    // truncate and normalize a value. Helper function for chargeDelta
    float truncNorm(float val, float max);
    
    // increment to a sphere's accumulated charge for a single timestep
    float chargeDelta(Eigen::Vector3f& loc, Eigen::Vector3f& vel, float density);
    
    // one time step increment of the charge for all spheres
    void update_sphere_charges();

    // Export
    void update_subspheres_params();
    // void export_spheres_combined(unsigned int plume_export_frequency);
    void export_ashfall(unsigned int plume_export_frequency);
    void export_spheres(unsigned int plume_export_frequency);
    void export_spheres_binary(unsigned int plume_export_frequency, bool headercount = false);

    // Fill structures
    void fill_height_field(std::vector<Eigen::Vector3f>& position, std::vector<Eigen::Vector3f>& normal,
                           Eigen::Matrix3f rotation, float scaling, Eigen::Vector3f translation);

    // Init
    void addAtmoModel(std::shared_ptr<PWM::Model::world<dsType, dsSType, valType, valType2>>& model);
    void setPlumeLoc(Eigen::Vector3f loc);
    //void set_gui();

    //gui_parameters gui_param;
};

#endif //CF_SMOKE_HPP
