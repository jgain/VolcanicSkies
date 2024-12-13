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
#ifndef PWM_SQUARE_2D_ARRAY_H
#define PWM_SQUARE_2D_ARRAY_H

#include "AbstractPWMDataStructure.h"
#include <functional>
#include "mathUtils.h"
#include <memory>
#include <omp.h>
#include <sstream>

namespace PWM{
    namespace PWMDataStructure{
        template<typename T> class square2DArray;
        
        template<typename T>
		void swap(square2DArray<T>& i, square2DArray<T>& j);

		template<typename T> class square2DArray: public AbstractPWMDataStructure<T>{
            private:
                size_t width = 0;
                double actualSize = 50000;
                Coordinate coord = Coordinate(0, 0);
            public:
                inline square2DArray();
                inline square2DArray(const size_t arrWidth);
                inline square2DArray(const size_t arrWidth, const double size);
                inline square2DArray(const square2DArray<T>& other);
                
                inline const T getData(const int i, const int j) const;
                inline const T getData(const size_t index) const;
                inline const T getInterpolated(const double x, const double y) const override;
                inline const size_t getWidth() const;
                inline const double getSize() const;
                
                inline const double gridLength() const override;
                inline const Coordinate getCoordinates() const;
                inline const Coordinate getCoordinates(const size_t index) const override;
                inline const Coordinate getCoordinates(const size_t i, const size_t j) const;

                inline const std::string print() const override;

                inline void setData(const size_t index, const T& val);
                inline void setData(const int i, const int j, const T& val);
                inline square2DArray& operator=(const T& val);
                inline void setCoord(const Coordinate& nCoord);

                inline void resize(const square2DArray<T>& target);
                inline void resize(const size_t arrWidth, const double size);
                inline void expand(square2DArray<T>& res, int factor = 2) const;
                inline void laplacian(square2DArray<T>& res) const;
                inline void gaussSeidelSmoothBS(square2DArray<T>& res, const square2DArray<T>& other, const T gridSize = 1.f);
                inline void gaussSeidelSmoothBSI(const square2DArray<T>& other, const int nb, const T gridSize = 1.f);

                inline bool operator==(const square2DArray<T>& other) const;
                inline bool operator!=(const square2DArray<T>& other) const;
                inline bool compareSizes(const square2DArray<T>& other) const;
        };

        template<typename T>
    inline square2DArray<T>::square2DArray() : AbstractPWMDataStructure<T>(), width(0){
        }

        template<>
    inline square2DArray<std::string>::square2DArray() : AbstractPWMDataStructure<std::string>(), width(0){
        }

        template<typename T>
    inline square2DArray<T>::square2DArray(const size_t arrWidth) : AbstractPWMDataStructure<T>(arrWidth * arrWidth), width(arrWidth){
            this->copy(0);
        }

        template<>
    inline square2DArray<std::string>::square2DArray(const size_t arrWidth) : AbstractPWMDataStructure<std::string>(arrWidth * arrWidth), width(arrWidth){
            this->copy("");
        }

        template<typename T>
    inline square2DArray<T>::square2DArray(const size_t arrWidth, const double size) : AbstractPWMDataStructure<T>(arrWidth * arrWidth), width(arrWidth), actualSize(size){
            this->copy(0);
        }

        template<>
    inline square2DArray<std::string>::square2DArray(const size_t arrWidth, const double size) : AbstractPWMDataStructure<std::string>(arrWidth * arrWidth), width(arrWidth), actualSize(size){
            this->copy("");
        }

        template<typename T>
    inline square2DArray<T>::square2DArray(const square2DArray<T>& other) : AbstractPWMDataStructure<T>(other.getWidth() * other.getWidth()), width(other.getWidth()), actualSize(other.getSize()){
            this->copy(other);
        }

        template<typename T>
    inline bool square2DArray<T>::operator==(const square2DArray<T>& other) const{
            if (getWidth() != other.getWidth())
                return false;
            if (getSize() != other.getSize())
                return false;
            for (int i = 0; i < getWidth(); ++i){
                for (int j = 0; j < getWidth(); ++j){
                    if (getData(i, j) != other.getData(i, j))
                        return false;
                }
            }
            return true;
        }

        template<typename T>
    inline bool square2DArray<T>::operator!=(const square2DArray<T>& other) const{
            return !(*this == other);
        }

        template<typename T>
    inline bool square2DArray<T>::compareSizes(const square2DArray<T>& other) const{
            if (getWidth() != other.getWidth())
                return false;
            if (getSize() != other.getSize())
                return false;
            return true;
        }

        template<typename T>
    inline const size_t square2DArray<T>::getWidth() const{
            return width;
        }

        template<typename T>
    inline const T square2DArray<T>::getData(const size_t index) const{
            return AbstractPWMDataStructure<T>::getData(index);
        }

        template<typename T>
    inline const T square2DArray<T>::getData(const int i, const int j) const{
            int i1 = i, j1 = j;
            if (i < 0)
                i1 = (int) (i % getWidth());
            else if (i >= getWidth())
                i1 = (int) (i % getWidth());
            if (j < 0)
                j1 = (int) (j % getWidth());
            else if (j >= getWidth())
                j1 = (int) (j % getWidth());
            return getData(i1 * getWidth() + j1);
        }

        template<typename T>
        inline const T square2DArray<T>::getInterpolated(const double x, const double y) const{
            float x1 = PWM::Utils::unsignedFloatingModulo(x, getWidth());
            float y1 = PWM::Utils::unsignedFloatingModulo(y, getWidth());
            float i, j, s, t;
            s = std::modf(x1, &i);
            t = std::modf(y1, &j);
			T a_ij = getData((int) i, (int) j);
			T a_i1j = getData((int) i + 1, (int) j);
			T a_ij1 = getData((int) i, (int) j + 1);
			T a_i1j1 = getData((int) i + 1, (int) j + 1);
			return PWM::Utils::interpolate(a_ij, a_i1j, a_ij1, a_i1j1, s, t);
        }

        template<>
        inline const std::string square2DArray<std::string>::getInterpolated(const double x, const double y) const{
            int i = std::round(x), j = std::round(y);
            return getData(i, j);
        }

        template<typename T>
    inline const double square2DArray<T>::getSize() const{
            return actualSize;
        }

        template<typename T>
    inline const double square2DArray<T>::gridLength() const{
            return actualSize / getWidth();
        }

        template<typename T>
    inline const Coordinate square2DArray<T>::getCoordinates() const{
            return coord;
        }

        template<typename T>
    inline const Coordinate square2DArray<T>::getCoordinates(const size_t index) const{
            return coord;
        }

        template<typename T>
    inline const Coordinate square2DArray<T>::getCoordinates(const size_t i, const size_t j) const{
            int i1 = i, j1 = j;
            if (i < 0)
                i1 = (int) (i % getWidth());
            else if (i >= getWidth())
                i1 = (int) (i % getWidth());
            if (j < 0)
                j1 = (int) (j % getWidth());
            else if (j >= getWidth())
                j1 = (int) (j % getWidth());
            return getCoordinates(i1 * getWidth() + j1);
        }

        template<typename T>
    inline const std::string square2DArray<T>::print() const{
            std::ostringstream r;
            r.precision(5);
            for (int i = 0; i < getWidth(); ++i){
                r << getData(i, 0);
                for (int j = 1; j < getWidth(); ++j){
                    r << '\t' << getData(i, j);
                }
                r << '\n';
            }
            return r.str();
        }

        template<typename T>
    inline void square2DArray<T>::setData(const size_t index, const T& val){
            AbstractPWMDataStructure<T>::setData(index, val);
        }

        template<typename T>
    inline void square2DArray<T>::setData(const int i, const int j, const T& val){
            int i1 = i, j1 = j;
            if (i < 0)
                i1 = (int) (i % getWidth());
            else if (i >= getWidth())
                i1 = (int) (i % getWidth());
            if (j < 0)
                j1 = (int) (j % getWidth());
            else if (j >= getWidth())
                j1 = (int) (j % getWidth());
            setData(i1 * getWidth() + j1, val);
        }

        template<typename T>
    inline square2DArray<T>& square2DArray<T>::operator=(const T& val){
            #pragma omp parallel for
                for (size_t i = 0; i < getWidth(); ++i)
                    for (size_t j = 0; j < getWidth(); ++j)
                        setData(i, j, val);
            #pragma omp barrier
            return *this;
        }

        template<typename T>
    inline void square2DArray<T>::setCoord(const Coordinate& nCoord){
            coord = nCoord;
        }

        template<typename T>
    inline void square2DArray<T>::resize(const square2DArray<T>& target){
            width = target.getWidth();
            actualSize = target.getSize();
            AbstractPWMDataStructure<T>::datasize = target.size();
        }

        template<typename T>
    inline void square2DArray<T>::expand(square2DArray<T>& res, int factor) const{
            if (res.getWidth() != getWidth() * factor){
                std::cerr << "Error in void square2DArray<T>::expand(square2DArray<T>& res, int factor) const: mismatched sizes for this and res.\nNo expansion performed." << std::endl;
                return;
            }
            #pragma omp parallel for
                for (int i = 0; i < getWidth() * factor; ++i){
                    for (int j = 0; j < getWidth() * factor; ++j)
                        res.setData(i, j, this->getData(i / factor, j / factor));
                }
            #pragma omp barrier
        }

        template<typename T>
    inline void square2DArray<T>::laplacian(square2DArray<T>& res) const{
            if (!compareSizes(res)){
                std::cerr << "Error in void square2DArray<T>::laplacian(square2DArray<T>& res) const: mismatched sizes for this and res.\nNo laplacian performed." << std::endl;
                return;
            }
            T h2 = 1.0 / (getWidth() * getWidth());

            #pragma omp parallel for
                for (int i = 0; i < getWidth(); ++i){
                    for (int j = 0; j < getWidth(); ++j){
                        T newVal = this->getData(i + 1, j) + this->getData(i - 1, j) + this->getData(i, j + 1) + this->getData(i, j - 1) - (4.0f * this->getData(i, j));
                        res.setData(i, j, newVal / h2);
                    }
                }
            #pragma omp barrier
        }

        template<typename T>
    inline void square2DArray<T>::gaussSeidelSmoothBS(square2DArray<T>& res, const square2DArray<T>& other, const T gridSize){
            if (getWidth() != res.getWidth() || getWidth() != other.getWidth()){
                std::cerr << "Error in void square2DArray<T>::gaussSeidelSmoothBS(square2DArray<T>& res, const square2DArray<T>& other, const T gridSize): mistmatched sizes for this and res or this and other.\nNo smoothing performed." << std::endl;
                return;
            }
            T h2 = (gridSize * gridSize) / (getWidth() * getWidth());

            #pragma omp parallel for
                for (int i = 0; i < getWidth(); ++i){
                    for (int j = 0; j < getWidth(); ++j){
                        T newVal = this->getData(i + 1, j) + this->getData(i - 1, j) + this->getData(i, j + 1) + this->getData(i, j - 1) - (h2 * other.getData(i, j));
                        res.setData(i, j, newVal / 4);
                    }
                }
            #pragma omp barrier
        }

        template<typename T>
    inline void square2DArray<T>::gaussSeidelSmoothBSI(const square2DArray<T>& other, const int nb, const T gridSize){
            if (getWidth() != other.getWidth()){
                std::cerr << "Error in void square2DArray<T>::gaussSeidelSmoothBS(const square2DArray<T>& other, const int nb, const T gridSize): mismatched sizes for this and other.\nNo smoothing performed." << std::endl;
                return;
            }
            T h2 = (gridSize * gridSize) / (getWidth() * getWidth());

            for (int i = 0; i < getWidth(); ++i){
                for (int j = 0; j < getWidth(); ++j){
                    T newVal = this->getData(i + 1, j) + this->getData(i - 1, j) + this->getData(i, j + 1) + this->getData(i, j - 1) - (h2 * other.getData(i, j));
                    setData(i, j, newVal / 4);
                }
            }
        }

        template<typename T>
    inline void swap(square2DArray<T>& i, square2DArray<T>& j){
			std::swap(i.faceWidth, j.faceWidth);
			std::swap((AbstractPWMDataStructure<T>) i, (AbstractPWMDataStructure<T>) j);
		}
    }
}
#endif //PWM_SQUARE_2D_ARRAY_H
