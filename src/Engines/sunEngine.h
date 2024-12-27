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
#ifndef PWM_SUN_ENGINE_H
#define PWM_SUN_ENGINE_H

#include "abstractEngine.h"
#include <cmath>
#include "sun.h"
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename V> class sunEngine : public AbstractEngine{
            private:
                std::vector<std::shared_ptr<PWM::Model::sun<V>>> sunList;

                void stepApparentRightAscension(std::shared_ptr<PWM::Model::sun<V>>& s);
                void stepApparentDeclination(std::shared_ptr<PWM::Model::sun<V>>& s);
            
            protected:
                void step_internal();
            
            public:
                sunEngine(float dt, bool active);
                const std::vector<std::shared_ptr<PWM::Model::sun<V>>>& getSunList() const;
                void addSun(std::shared_ptr<PWM::Model::sun<V>>& s);
        };

        template<typename V>
        inline sunEngine<V>::sunEngine(float dt, bool active) : AbstractEngine(dt, active){
            sunList = std::vector<std::shared_ptr<PWM::Model::sun<V>>>();
        }

        template<typename V>
        inline void sunEngine<V>::addSun(std::shared_ptr<PWM::Model::sun<V>>& s){
            sunList.push_back(s);
        }

        template<typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::sun<V>>>& sunEngine<V>::getSunList() const{
            return sunList;
        }

        template<typename V>
        inline void sunEngine<V>::step_internal(){
            for (auto s : sunList){
                stepApparentDeclination(s);
                stepApparentRightAscension(s);
            }
        }

        template<typename V>
        inline void sunEngine<V>::stepApparentRightAscension(std::shared_ptr<PWM::Model::sun<V>>& s){
            if (s->getApparentRightAscension() > 180.){
                V x = s->getApparentRightAscension() - 360;
                s->setApparentRightAscension(x);
            }
            V r = s->getApparentRightAscension() + (this->getDt() * (360 / s->getSolarDayLength()));
            s->setApparentRightAscension(r);
        }

        template<typename V>
        inline void sunEngine<V>::stepApparentDeclination(std::shared_ptr<PWM::Model::sun<V>>& s){
            if (std::abs(s->getApparentDeclination()) > s->getAxialTilt())
                s->flipSeasonDirection();
            V r = s->getApparentDeclination() + (this->getDt() * ((4 * s->getAxialTilt()) / s->getTropicalYearLength()) * s->getSeasonDirection());
            s->setApparentDeclination(r);
        }
    }
}

#endif //PWM_SUN_ENGINE_H
