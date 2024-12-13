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
