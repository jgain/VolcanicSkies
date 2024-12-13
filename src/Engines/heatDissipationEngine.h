#ifndef PWM_HEAT_DISSIPATION_ENGINE_H
#define PWM_HEAT_DISSIPATION_ENGINE_H

#include "abstractEngine.h"
#include "airLayer.h"
#include <cmath>
#include <omp.h>
#include <vector>

namespace PWM{
    namespace Engine{
        template<typename T, typename V> class heatDissipationEngine : public AbstractEngine{
            private:
                V coefficient;
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;
            protected:
                void step_internal();
            public:
                heatDissipationEngine(float dt, bool active, V coeff = -1);

                void addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l);
                const V& getCoefficient() const;
                void setCoefficient(V val);
        };

        template<typename T, typename V>
        inline heatDissipationEngine<T, V>::heatDissipationEngine(float dt, bool active, V coeff) : AbstractEngine(dt, active), coefficient(coeff){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            if (coefficient == -1)
                coefficient = 0.001 * dt;
        }

        template<typename T, typename V>
        inline void heatDissipationEngine<T, V>::addAirLayer(std::shared_ptr<PWM::Model::airLayer<T, V>> l){
            airLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const V& heatDissipationEngine<T, V>::getCoefficient() const{
            return coefficient;
        }

        template<typename T, typename V>
        inline void heatDissipationEngine<T, V>::setCoefficient(V val){
            coefficient = val;
        }

        template<typename T, typename V>
        inline void heatDissipationEngine<T, V>::step_internal(){
            V x = std::exp(-getDt() * coefficient);
            for (auto l : airLayers)
                l->getTemperature().mul(x);
        }
    }
}
#endif //PWM_HEAT_DISSIPATION_ENGINE_H
