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
#ifndef PWM_MODEL_CONVECTION_LAYER_H
#define PWM_MODEL_CONVECTION_LAYER_H

#include <iostream>
#include "planet.h"

namespace PWM{
    namespace Model{
        template<typename T, typename V> class convectionLayer{
            private:
                //The height at the top of this convectionLayer in metres.
                V topHeight;

                //The height at the bottom of this convectionLayer in metres.
                V bottomHeight;

                V height;
                //The vertical velocities in a grid structure, in metres per second
                std::shared_ptr<T> verticalVelocities;

                //The rainfall in this layer in a grid structure, in kg per cubic metre
                std::shared_ptr<T> rainfall;

                //The ashfall in this layer in a grid structure, in kg per cubic metre
                std::shared_ptr<T> ashfall;

                std::shared_ptr<planet> Planet;
            public:
                //Default constructor not meant to be used
                convectionLayer(std::shared_ptr<planet> p, const V bot, const V top);

                //Constructor meant for use with square2DArray
                convectionLayer(std::shared_ptr<planet> p, const V bot, const V top, const size_t faceWidth, const V size = 50000.f);

                //Constructor meant for use with square2DArray
                convectionLayer(std::shared_ptr<planet> p, const V bot, const V top, const size_t xWidth, const size_t yHeight, const V xSize, const V ySize);

                const std::shared_ptr<planet> getPlanet() const;
                
                const T& getVerticalVelocities() const;
                const T& getRainfall() const;
                const T& getAshfall() const;

                const V getVerticalVelocity(const size_t index) const;
                const V getVerticalVelocity(const size_t i, const size_t j) const;
                const V getRainfall(const size_t index) const;
                const V getRainfall(const size_t i, const size_t j) const;
                const V getAshfall(const size_t index) const;
                const V getAshfall(const size_t i, const size_t j) const;

                const V getBottomHeight() const;
                const V getTopHeight() const;
                const V getLayerHeight() const;
                const V getThickness() const;

                void setVerticalVelocity(const size_t index, const V& newVal);
                void setVerticalVelocity(const size_t i, const size_t j, const V& newVal);
                void setRainfall(const size_t index, const V& newVal);
                void setRainfall(const size_t i, const size_t j, const V& newVal);
                void setAshfall(const size_t index, const V& newVal);
                void setAshfall(const size_t i, const size_t j, const V& newVal);

                void swapVertVel(std::shared_ptr<T>& oldVertVels);

                bool operator==(const convectionLayer<T, V>& other) const;
                bool operator!=(const convectionLayer<T, V>& other) const;

                template<typename X, typename Y>
                friend bool operator<(const convectionLayer<X, Y>& lhs, const convectionLayer<X, Y>& rhs);
                template<typename X, typename Y>
                friend bool operator<(std::shared_ptr<convectionLayer<X, Y>> lhs, std::shared_ptr<convectionLayer<X, Y>> rhs);
        };

        template<typename T, typename V>
    inline convectionLayer<T, V>::convectionLayer(std::shared_ptr<planet> p, const V bot, const V top) : Planet(p), topHeight(top), bottomHeight(bot){
            verticalVelocities = nullptr;
            rainfall = nullptr;
            ashfall = nullptr;
            height = static_cast<V>((top + bot) / 2);
        }

        template<typename T, typename V>
    inline convectionLayer<T, V>::convectionLayer(std::shared_ptr<planet> p, const V bot, const V top, const size_t faceWidth, const V size) : Planet(p), bottomHeight(bot), topHeight(top){
            verticalVelocities = std::make_shared<T>(faceWidth, size);
            rainfall = std::make_shared<T>(faceWidth, size);
            ashfall = std::make_shared<T>(faceWidth, size);
            verticalVelocities->copy(0.);
            rainfall->copy(0.);
            ashfall->copy(0);
            height = static_cast<V>((top + bot) / 2);
        }

        template<typename T, typename V>
    inline convectionLayer<T, V>::convectionLayer(std::shared_ptr<planet> p, const V bot, const V top, const size_t xWidth, const size_t yHeight, const V xSize, const V ySize) : Planet(p), bottomHeight(bot), topHeight(top){
            verticalVelocities = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            rainfall = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            ashfall = std::make_shared<T>(xWidth, yHeight, 0.5, 0.5, xSize, ySize);
            verticalVelocities->copy(0.);
            rainfall->copy(0.);
            ashfall->copy(0.);
            height = static_cast<V>((top + bot) / 2);
        }

        template<typename T, typename V>
    inline const std::shared_ptr<planet> convectionLayer<T, V>::getPlanet() const{
            return Planet;
        }

        template<typename T, typename V>
    inline const T& convectionLayer<T, V>::getVerticalVelocities() const{
            return *verticalVelocities;
        }

        template<typename T, typename V>
    inline const T& convectionLayer<T, V>::getRainfall() const{
            return *rainfall;
        }

        template<typename T, typename V>
    inline const T &convectionLayer<T, V>::getAshfall() const{
            return *ashfall;
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getVerticalVelocity(const size_t index) const{
            return verticalVelocities->getData(index);
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getVerticalVelocity(const size_t i, const size_t j) const{
            return verticalVelocities->getData(i, j);
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getRainfall(const size_t index) const{
            return rainfall->getData(index);
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getRainfall(const size_t i, const size_t j) const{
            return rainfall->getData(i, j);
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getAshfall(const size_t index) const{
            return ashfall->getData(index);
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getAshfall(const size_t i, const size_t j) const{
            return ashfall->getData(i, j);
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getBottomHeight() const{
            return bottomHeight;
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getTopHeight() const{
            return topHeight;
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getLayerHeight() const{
            return height;
        }

        template<typename T, typename V>
    inline const V convectionLayer<T, V>::getThickness() const{
            return topHeight - bottomHeight;
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::setVerticalVelocity(const size_t index, const V& newVal){
            verticalVelocities->setData(index, newVal);
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::setVerticalVelocity(const size_t i, const size_t j, const V& newVal){
            verticalVelocities->setData(i, j, newVal);
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::setRainfall(const size_t index, const V& newVal){
            rainfall->setData(index, newVal);
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::setRainfall(const size_t i, const size_t j, const V& newVal){
            rainfall->setData(i, j, newVal);
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::setAshfall(const size_t index, const V& newVal){
            ashfall->setData(index, newVal);
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::setAshfall(const size_t i, const size_t j, const V& newVal){
            ashfall->setData(i, j, newVal);
        }

        template<typename T, typename V>
    inline void convectionLayer<T, V>::swapVertVel(std::shared_ptr<T> &oldVertVels){
            std::swap(verticalVelocities, oldVertVels);
        }

        template<typename T, typename V>
    inline bool convectionLayer<T, V>::operator==(const convectionLayer<T, V>& other) const{
            if (getPlanet() != other.getPlanet())
                return false;
            if (getTopHeight() != other.getTopHeight())
                return false;
            if (getBottomHeight() != other.getBottomHeight())
                return false;
            if (getVerticalVelocities() != other.getVerticalVelocities())
                return false;
            if (getRainfall() != other.getRainfall())
                return false;
            return true;
        }

        template<typename T, typename V>
    inline bool convectionLayer<T, V>::operator!=(const convectionLayer<T, V>& other) const{
            return !(*this == other);
        }

        template<typename T, typename V>
    inline bool operator<(const convectionLayer<T, V>& lhs, const convectionLayer<T, V>& rhs){
            return lhs.getLayerHeight() < rhs.getLayerHeight();
        }

        template<typename T, typename V>
    inline bool operator<(std::shared_ptr<convectionLayer<T, V>> lhs, std::shared_ptr<convectionLayer<T, V>> rhs){
            return lhs->getLayerHeight() < rhs->getLayerHeight();
        }
    }
}

#endif //PWM_MODEL_CONVECTION_LAYER_H
