#ifndef CF_SMOKE_LAYER_HPP
#define CF_SMOKE_LAYER_HPP

#include <Eigen/Dense>
#include <Eigen/Geometry>



struct subsphere_params
{
    size_t parent_id;
    Eigen::Vector3f relative_position;

    Eigen::Vector3f center;
    Eigen::Vector3f speed;
    float r;
    float charge; // for volcanic lightning
    float size_ratio;
    float z_factor; // for ellipsoids
    float minz, maxz; // for skirt intersection

    subsphere_params() {}
};

struct free_sphere_params
{
    // characteristics
    Eigen::Vector3f center;
    float r;

    Eigen::Vector3f speed;
    float mass;
    float rho;
    float charge; // for volcanic lightning
    float relative_distance;
    float minz, maxz; // for skirt intersection

    float density_loss_last_dt;

    // rotation informations
    float angular_speed;
    Eigen::Vector3f rotation_axis; // cst
    float current_angle; // to keep track of self rotation
    Eigen::Vector3f angle_vector;

    // for diversity
    float perturbation;
    float size_factor;

    // ellipsoids in stagnating spheres
    float z_factor;

    // identify state of spheres
    bool falling_under_atm_rho;
    bool stagnate;
    bool stagnate_long;
    bool falling;
    bool falling_disappeared; // to make their disappearance more discrete
    bool secondary_column;

    // informations for stagnation
    int closest_layer_idx;
    float stagnation_altitude;
    float max_altitude;
    float xy_at_max_altitude;
    Eigen::Vector3f center_at_max_altitude;
    Eigen::Vector3f layer_center_at_max_altitude;


    free_sphere_params() {}
    // constructor used for spheres attached to newly emitted slice
    free_sphere_params(Eigen::Vector3f ring_center, float angle_on_ring, float r, float speed, float charge, bool sec_column=false) : r(r), speed({0,0,speed}), charge(charge), secondary_column(sec_column)
    {
        Eigen::Vector3f angle_normal = {cos(angle_on_ring), sin(angle_on_ring),0};
        angle_vector = angle_normal;
        center = ring_center + r*angle_normal;
        rotation_axis = {-sin(angle_on_ring), cos(angle_on_ring),0};
        angular_speed = speed/r;
        current_angle = 0;
        perturbation = 0;
        z_factor = 1.f;
        stagnate = false;
        stagnate_long = false;
        falling = false;
        falling_disappeared = false;
        density_loss_last_dt = 0.;
    }
    // constructor used for falling spheres
    free_sphere_params(Eigen::Vector3f center, float r, float rho, float charge): r(r), center(center), rho(rho), charge(charge), speed({0,0,0}), angular_speed(0), rotation_axis({1,0,0}), current_angle(0),
        angle_vector({0,1,0}), mass(rho*4./3.*3.14*r*r*r), size_factor(1.), falling_under_atm_rho(false), stagnate(false), stagnate_long(false), falling(true), falling_disappeared(false) {}
};

struct smoke_layer
{
    Eigen::Vector3f center; // Position (x,y,z) (m)
    Eigen::Vector3f v; // Speed (m.s-1)
    Eigen::Vector3f a; // Acceleration (m.s-2)

    Eigen::Vector3f plume_axis; // actually normalized speed
    float speed_along_axis; // (m.s-1)
    float theta; // angle of plume axis (rad)
    Eigen::Vector3f theta_axis; // axis to compute plume_axis rotation with theta

    float r; // ray (m)
    float temperature; // (Kelvin)
    float moisture; // absolute humidity (kg.m-3)
    float rho; // density (kg.m-3)
    float charge; // static electricity proxy, maximum charge per time-step = 1
    float thickness; // (m)
    float minz, maxz; // for skirt intersection

    // identify state of layers
    bool rising;
    bool begin_falling; // true for first frame of falling to allow closest spheres to fall
    bool falling;
    bool plume; // becomes true if rho becomes lower than atm
    bool stagnates; // becomes true if densities of layer and air become equal
    bool stagnates_long; // becomes true if stagnates and stops rising
    bool secondary_plume;


    smoke_layer() : center{0,0,0}, v{0,0,200}, a{0,0,0}, speed_along_axis(200), theta(asin(1.0)), plume_axis{0,0,1}, theta_axis{1,0,0},
        r(100), temperature(1273), rho(1.5), charge(3.0), thickness(100),
        rising(true), begin_falling(false), falling(false), plume(false), stagnates(false), stagnates_long(false) {}
    // constructor for principal plume layers
    smoke_layer(Eigen::Vector3f v) : center{0,0,0}, v(v), a{0,0,0}, speed_along_axis(v.norm()), theta(asin(1.0)), plume_axis(v.normalized()),
        r(100), temperature(1273), rho(1.5), charge(3.0), thickness(100),
        rising(true), begin_falling(false), falling(false), plume(false), stagnates(false), stagnates_long(false) {
        if (v.norm() == v.z()) theta_axis = Eigen::Vector3f(1,0,0);
        else theta_axis = (v.cross(Eigen::Vector3f(0,0,1))).normalized();
    }
    smoke_layer(Eigen::Vector3f p, Eigen::Vector3f v) : center(p), v(v), a{0,0,0}, speed_along_axis(v.norm()), theta(asin(1.0)), plume_axis(v.normalized()),
        r(100), temperature(1273), rho(1.5), charge(3.0), thickness(100),
        rising(true), begin_falling(false), falling(false), plume(false), stagnates(false), stagnates_long(false) {
        if (v.norm() == v.z()) theta_axis = Eigen::Vector3f(1,0,0);
        else theta_axis = (v.cross(Eigen::Vector3f(0,0,1))).normalized();
    }
    // constructor for secondary plume layers
    smoke_layer(Eigen::Vector3f v, float rho, float chg, float r, Eigen::Vector3f position, bool secondary_plume) : center(position), v(v), a{0,0,0}, speed_along_axis(v.norm()), theta(asin(1.0)), plume_axis(v.normalized()),
        r(r), temperature(1273), rho(rho), charge(chg), thickness(r),
        rising(true), begin_falling(false), falling(false), plume(false), stagnates(false), stagnates_long(false), secondary_plume(secondary_plume) {
        if (v.norm() == v.z()) theta_axis = Eigen::Vector3f(1,0,0);
        else theta_axis = (v.cross(Eigen::Vector3f(0,0,1))).normalized();
    }
};

#endif //CF_SMOKE_LAYER_HPP
