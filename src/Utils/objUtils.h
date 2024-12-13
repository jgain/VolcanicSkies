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
#ifndef OBJ_UTILS_H
#define OBJ_UTILS_H

#include <Eigen/Dense>
#include "OBJ_Loader.h"

namespace PWM{
    namespace Utils{
    inline std::vector<Eigen::Vector3f> getPositions(std::string file){
            objl::Loader obj;
            obj.LoadFile(file);
            std::vector<Eigen::Vector3f> res;
            for (auto i : obj.LoadedVertices){
                Eigen::Vector3f point(i.Position.X, i.Position.Y, i.Position.Z);
                res.push_back(point);
            }
            return res;
        }

    inline std::vector<Eigen::Vector3f> getNormals(std::string file){
            objl::Loader obj;
            obj.LoadFile(file);
            std::vector<Eigen::Vector3f> res;
            for (auto i : obj.LoadedVertices){
                Eigen::Vector3f point(i.Normal.X, i.Normal.Y, i.Normal.Z);
                res.push_back(point);
            }
            return res;
        }

    inline std::pair<std::vector<Eigen::Vector3f>, std::vector<Eigen::Vector3f>> getPosAndNorms(std::string file){
            objl::Loader obj;
            obj.LoadFile(file);
            std::vector<Eigen::Vector3f> pos, norm;
            for (auto i : obj.LoadedVertices){
                Eigen::Vector3f posI(i.Position.X, i.Position.Y, i.Position.Z);
                pos.push_back(posI);
                Eigen::Vector3f normI(i.Normal.X, i.Normal.Y, i.Normal.Z);
                norm.push_back(normI);
            }
            return std::make_pair(pos, norm);
        }
    }
}
#endif //OBJ_UTILS_H
