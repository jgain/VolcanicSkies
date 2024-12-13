#ifndef CF_TERRAIN_HPP
#define CF_TERRAIN_HPP

#include <Eigen/Dense>
#include <vector>

// Terrain grid for acceleration collision computation
struct terrain_structure
{
    std::vector<Eigen::Vector3f> positions;
    std::vector<Eigen::Vector3f> normals;

    std::vector< std::vector< std::vector<unsigned int> > > grid;
    size_t grid_size;

    std::vector< std::vector<float> > height_field;
    //vcl::mesh_drawable height_field_mesh;
    std::vector< std::vector<Eigen::Vector3f> > normal_field;
    float cell_size;
    size_t field_size;
    float min_xyz;
    float max_xyz;
    float latitude;

    terrain_structure();
    terrain_structure(std::string file);
    terrain_structure(std::string file, float cellSize, float sceneWidth);
    void constructTerrainFromObjFile(std::string file);
    void constructTerrainFromElvFile(std::string file);
    terrain_structure& testTerrainStructure(int width);
    void fill_height_field(std::vector<Eigen::Vector3f>& positions, std::vector<Eigen::Vector3f>& normals, Eigen::Matrix3f rotation, float scaling, Eigen::Vector3f translation);
    void write_height_field(std::string file);
    float sampleHeightAt(float i, float j);
};

#endif //CF_TERRAIN_HPP
