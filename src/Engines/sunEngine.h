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
