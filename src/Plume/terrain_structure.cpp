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
#include <cmath>
#include "mathUtils.h"
#include "objUtils.h"
#include "terrain_structure.hpp"

terrain_structure::terrain_structure(){
    height_field = std::vector<std::vector<float>>();
}

terrain_structure::terrain_structure(std::string file){
    this->constructTerrainFromElvFile(file);
}

terrain_structure::terrain_structure(std::string file, float cellSize, float sceneWidth){
    // grid_size = someNumber?
    this->min_xyz = -sceneWidth;
    this->max_xyz = sceneWidth;
    this->cell_size = cellSize;
    this->field_size = (size_t)((2 * this->max_xyz) / this->cell_size);
    this->latitude = 0.0;
    this->constructTerrainFromObjFile(file);
}

void terrain_structure::constructTerrainFromObjFile(std::string file){
    auto v = PWM::Utils::getPosAndNorms(file);
    this->positions = v.first;
    this->normals = v.second;

    Eigen::Matrix3f rotation = (Eigen::AngleAxisf(3.14f/2.0f, Eigen::Vector3f(1.0f, 0, 0))).toRotationMatrix();
    float scaling = 100.f;
    Eigen::Vector3f translation(37.f,68.f,-25.f);
    this->fill_height_field(positions, normals, rotation, scaling, translation);
    write_height_field("/Users/james/Desktop/mtsthelens.elv");
}

void terrain_structure::constructTerrainFromElvFile(std::string file){
    try {
        auto fil = std::ifstream(file);
        int width, height;
        float cellSize, lat;

        //Read header
        fil >> width;
        fil >> height;
        fil >> cellSize;
        fil >> lat;

        this->field_size = width;
        this->cell_size = cellSize;
        this->latitude = lat;

        this->min_xyz = -cellSize * (field_size * 0.5);
        this->max_xyz = -this->min_xyz;

        height_field.clear();
        //Read grid points
        for (int i = 0; i < height; ++i){
            height_field.push_back(std::vector<float>());
            for (int j = 0; j < width; ++j){
                float val;
                fil >> val;
                height_field[i].push_back(val);
            }
        }
        normal_field.clear();
        for (int i = 0; i < height; ++i){
            normal_field.push_back(std::vector<Eigen::Vector3f>());
            for (int j = 0; j < width; ++j){
                float aB, aT, aR, aL;
                if (i == 0)
                    aT = height_field[i][j];
                else
                    aT = height_field[i - 1][j];
                if (i == height - 1)
                    aB = height_field[i][j];
                else
                    aB = height_field[i + 1][j];
                if (j == 0)
                    aL = height_field[i][j];
                else
                    aL = height_field[i][j - 1];
                if (j == width - 1)
                    aR = height_field[i][j];
                else
                    aR = height_field[i][j + 1];
                Eigen::Vector3f norm = Eigen::Vector3f((aR - aL) / (2 * this->cell_size), (aB - aT) / (2 * this->cell_size), -1).normalized();
                normal_field[i].push_back(norm);
            }
        }
    } catch (std::exception& e) {
        std::cerr << "Error trying to read file " << file << " in as elv terrain." << std::endl;
        assert(0);
    }

}

void terrain_structure::write_height_field(std::string file){
    std::ofstream outfile;

    outfile.open((char *) file.c_str(), std::ios_base::out);
    if(outfile.is_open())
    {
        outfile << this->field_size << " " << this->field_size << " " << this->cell_size << " " << 46.19 << std::endl; // latitude of mt st helens = 46.19
        for (unsigned int i = 0; i<this->field_size; i++)
        {
            for (unsigned int j = 0; j<this->field_size; j++)
            {
                outfile << this->height_field[i][j] / 0.3048f << " ";
            }
        }
        outfile << std::endl;
        outfile.close();
    }
    else
    {
        std::cerr << "Error Terrain::loadElv:unable to open file " << file << std::endl;
    }
}

float terrain_structure::sampleHeightAt(float i, float j){
    float x = i / cell_size;
    float y = j / cell_size;
    int xIndex = (int) std::floor(x);
    int yIndex = (int) std::floor(y);
    double alphaX = x - xIndex, alphaY = y - yIndex;
    int x1, x2, y1, y2;
    x1 = xIndex % field_size;
    x2 = (x1 + 1) % field_size;
    y1 = yIndex % field_size;
    y2 = (y1 + 1) % field_size;

    float val11, val12, val21, val22;
    val11 = height_field[x1][y1];
    val12 = height_field[x1][y2];
    val21 = height_field[x2][y1];
    val22 = height_field[x2][y2];
    return PWM::Utils::interpolate(val11, val12, val21, val22, alphaY, alphaX);
}

terrain_structure& terrain_structure::testTerrainStructure(int width){
    float sceneWidth = 15000;
    this->field_size = width;
    this->min_xyz = -sceneWidth;
    this->max_xyz = sceneWidth;
    this->cell_size = 2 * sceneWidth / this->field_size;

    this->height_field.resize(this->field_size);
    this->normal_field.resize(this->field_size);
    for (unsigned int i = 0; i < this->field_size; ++i){
        this->height_field[i].resize(this->field_size);
        this->normal_field[i].resize(this->field_size);
    }
//    Eigen::Matrix3f rotation = (Eigen::AngleAxisf(3.14f/2.0f, Eigen::Vector3f(0, 0, 1.0f))).toRotationMatrix();
    for (int i = 0; i < this->field_size; ++i){
        for (int j = 0; j < this->field_size; ++j){
            this->height_field[i][j] = 0;
//            this->normal_field[i][j] = this->normal_field[i][j] * rotation;
        }
    }
    return *this;
}

void terrain_structure::fill_height_field(std::vector<Eigen::Vector3f>& position, std::vector<Eigen::Vector3f>& normal,
                                    Eigen::Matrix3f rotation, float scaling, Eigen::Vector3f translation)
{
    // prepare field with parameters
    float interval_size = this->max_xyz - this->min_xyz;
    this->height_field.resize(this->field_size);
    this->normal_field.resize(this->field_size);
    for (unsigned int i = 0; i < this->field_size; i++){
        this->height_field[i].resize(this->field_size);
        this->normal_field[i].resize(this->field_size);
    }

    //transform like for mesh_drawable
    for (unsigned int i = 0; i < position.size(); i++){
        position[i] = scaling * (rotation * position[i] + translation);
        normal[i] = rotation * normal[i];
    }
    this->positions = position;
    this->normals = normal;

    // test: check min and max x and y for mesh -> min_xyz should be chosen bigger than min_x and min_y
    /*float min_x = 0, min_y = 0, max_x = 0, max_y = 0;
    for (unsigned int i = 0; i < position.size(); i++)
    {
        if (position[i].x() < min_x)
        {
            min_x = position[i].x();
        }
        if (position[i].x() > max_x)
        {
            max_x = position[i].x();
        }
        if (position[i].y() < min_y)
        {
            min_y = position[i].y();
        }
        if (position[i].y() > max_y)
        {
            max_y = position[i].y();
        }
    }
    std::cout << "min x " << min_x << std::endl;
    std::cout << "max x " << max_x << std::endl;
    std::cout << "min y " << min_y << std::endl;
    std::cout << "max y " << max_y << std::endl;*/

    // find minimum height (for water regions like lakes)
    float min_height = position[0].z();
    unsigned int min_idx = 0;
    for (unsigned int i = 0; i < position.size(); i++)
    {
        if (position[i].z() < min_height)
        {
            min_idx = i;
            min_height = position[i].z();
        }
    }
    // initialize height field with this minimal value
    for (unsigned int i = 0; i<this->field_size; i++)
    {
        for (unsigned int j = 0; j<this->field_size; j++)
        {
            this->height_field[i][j] = min_height;
            this->normal_field[i][j] = Eigen::Vector3f(0,0,1);
        }
    }

    // pre computation for grid filling
    unsigned int mesh_pts_nb = position.size();
    unsigned int mesh_pts_in_one_row = sqrt((float)mesh_pts_nb); // works only if the terrain is a square
    unsigned int grid_cells_in_one_row = this->field_size;
    float nb_cells_between_two_pts = grid_cells_in_one_row / mesh_pts_in_one_row;

    // fill height field for collisions
    std::vector< std::vector<bool> > filled;
    filled.resize(this->field_size);
    for (unsigned int i = 0; i<this->field_size; i++)
    {
        filled[i].resize(this->field_size);
        for (unsigned int j = 0; j<this->field_size; j++)
        {
            filled[i][j] = false;
        }
    }
    for (unsigned int i = 0; i < position.size(); i++)
    {
        if (position[i].x() < this->max_xyz && position[i].x() > this->min_xyz && position[i].y() < this->max_xyz && position[i].y() > this->min_xyz)
        {
            int idx_x = (int)(this->field_size * (position[i].x() - this->min_xyz) / interval_size);
            int idx_y = (int)(this->field_size * (position[i].y() - this->min_xyz) / interval_size);
            if (idx_x == this->field_size) idx_x = this->field_size - 1;
            if (idx_y == this->field_size) idx_y = this->field_size - 1;
            this->height_field[idx_x][idx_y] = position[i].z();
            this->normal_field[idx_x][idx_y] = normal[i];
            filled[idx_x][idx_y] = true;

            if (nb_cells_between_two_pts > 1)
                for (unsigned int j = 0; j < nb_cells_between_two_pts + 1; j++)
                    for (unsigned int k = 0; k < nb_cells_between_two_pts + 1; k++)
                    {
                        this->height_field[idx_x+j][idx_y+k] = position[i].z();
                        this->normal_field[idx_x+j][idx_y+k] = normal[i];
                    }
        }
    }

    // test code to see which cells are not filled
    /*std::cout << "unfilled cells" << std::endl;
    unsigned int unfilled_cells = 0;
    for (unsigned int i = 0; i<this->field_size; i++)
    {
        for (unsigned int j = 0; j<this->field_size; j++)
        {
            if (filled[i][j] == false)
            {
                std::cout << "0 ";
            }
            else
                std::cout << "1 ";
        }
        std::cout << std::endl;
    }
    for (unsigned int i = 0; i<this->field_size; i++)
    {
        for (unsigned int j = 0; j<this->field_size; j++)
        {
            if (filled[i][j] == false)
            {
                std::cout << i << " " << j << std::endl;
                unfilled_cells++;
            }
        }
    }
    std::cout << "unfilled cells " << unfilled_cells << std::endl;*/
}
