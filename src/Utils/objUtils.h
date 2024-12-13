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
