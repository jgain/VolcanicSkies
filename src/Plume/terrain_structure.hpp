/*
 * Volcanic Skies
 * Copyright (C) 2024 P. Cilliers Pretorius, University of Cape Town
 *
 * This file is part of the Volcanic Skies project.
 *
 * Volcanic Skies is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License (GPL) as published by the Free Software 
 * Foundation, either version 2 of the License, or (at your discretion) any later version.
 *
 * Volcanic Skies is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with 
 * Volcanic Skies in the LICENSE file. If not, see <https://www.gnu.org/licenses/>.
 *
 * Additional information and disclaimers regarding liability and third-party 
 * components can be found in the NOTICE file included with this project.
 *
 */
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
