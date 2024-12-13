#ifndef CUBED_SPHERE_DS_H
#define CUBED_SPHERE_DS_H

#include "AbstractPWMDataStructure.h"
#include <cmath>
#include "mathUtils.h"

namespace PWM{
    namespace PWMDataStructure{
        template <typename T> class cubedSphereDS;

		template<typename T>
		void swap(cubedSphereDS<T>& i, cubedSphereDS<T>& j);

		template<typename T>
		class cubedSphereDS : public AbstractPWMDataStructure<T>{
			private:
				//Each face of the cube is N * N cells 
				size_t faceWidth;
			
			public:
				cubedSphereDS(size_t n = 4);
				cubedSphereDS(const cubedSphereDS& other);
				/*cubedSphereDS& operator=(const cubedSphereDS& other);*/

				bool operator==(const cubedSphereDS& other) const;
				bool operator!=(const cubedSphereDS& other) const;

				cubedSphereDS& operator=(const T& val); 

				const T getData(const size_t index) const;
				const T getData(const size_t face, const size_t i, const size_t j) const;
				
				void setData(const size_t index, const T& val);
				void setData(const size_t face, const size_t i, const size_t j, const T& val);
				
				const Coordinate getCoordinates(const size_t index) const override;
				const Coordinate getCoordinates(const size_t face, const size_t i, const size_t j) const;
				const Coordinate getEquatorialFaceCoords(const double leftEdgeLongitude, const size_t i, const size_t j) const;
				const Coordinate getPolarFaceCoords(const bool NP, const size_t i, const size_t j) const;

				const size_t getFaceWidth() const;
				const T gridLength() const override;

				const std::tuple<size_t, size_t, size_t> convertIndexToFaceCoord(const size_t index) const;
				const size_t convertFaceCoordToIndex(const size_t face, const size_t i, const size_t j) const;

				friend void swap <> (cubedSphereDS<T>& i, cubedSphereDS<T>& j);
		};

		template<typename T>
		cubedSphereDS<T>::cubedSphereDS(size_t n): AbstractPWMDataStructure<T>(6 * n * n){
			this->faceWidth = n;
		}

		template<typename T>
		cubedSphereDS<T>::cubedSphereDS(const cubedSphereDS<T>& other): AbstractPWMDataStructure<T>(other.size()){
			this->faceWidth = other.getFaceWidth();
			this->copy(other);
		}

		/*template<typename T>
		cubedSphereDS<T>& cubedSphereDS<T>::operator=(const cubedSphereDS<T>& other){
			*this()
			return *this;
		}*/

		template<typename T>
		bool cubedSphereDS<T>::operator==(const cubedSphereDS<T> & other) const {
			if (typeid(getData(0)) != typeid(other.getData(0))){
				return false;
			}
			if (!this->checkSize(other)){
				return false;
			}
			for (int i = 0; i < this->size(); ++i){
				if (getData(i) != other.getData(i)) {
					return false;
				}
			}
			return true;
		}

		template<typename T>
		bool cubedSphereDS<T>::operator!=(const cubedSphereDS<T> & other) const {
			return !(*this == other);
		}

		template<typename T>
		cubedSphereDS<T>& cubedSphereDS<T>::operator=(const T& val){
			AbstractPWMDataStructure<T>::operator=(val);
			return *this;
		}

		template<typename T>
		const T cubedSphereDS<T>::getData(const size_t index) const{
			return AbstractPWMDataStructure<T>::getData(index);
		}

		template<typename T>
		const T cubedSphereDS<T>::getData(const size_t face, const size_t i, const size_t j) const{
			return getData(convertFaceCoordToIndex(face, i, j));
		}

		template<typename T>
		void cubedSphereDS<T>::setData(const size_t index, const T& val){
			return AbstractPWMDataStructure<T>::setData(index, val);
		}

		template<typename T>
		void cubedSphereDS<T>::setData(const size_t face, const size_t i, const size_t j, const T& val){
			setData(convertFaceCoordToIndex(face, i, j), val);
		}

		template<typename T>
		const size_t cubedSphereDS<T>::getFaceWidth() const {
			return faceWidth;
		}

		template<typename T>
		const T cubedSphereDS<T>::gridLength() const {
			#if __cplusplus > 201703L
				double pi = std::numbers::pi;
			#else
				double pi = 3.1415926535;
			#endif
			return pi / (2 * getFaceWidth());
		}

		template<typename T>
		const Coordinate cubedSphereDS<T>::getCoordinates(const size_t index) const{
			auto a = convertIndexToFaceCoord(index);
			return getCoordinates(std::get<0>(a), std::get<1>(a), std::get<2>(a));
		}

		template<typename T>
		const Coordinate cubedSphereDS<T>::getCoordinates(const size_t face, const size_t i, const size_t j) const{
			switch(face){
				case 0:
					return getPolarFaceCoords(true, i, j);
					break;
				case 1:
					return getEquatorialFaceCoords(-180, i, j);
					break;
				case 2:
					return getEquatorialFaceCoords(-90, i, j);
					break;
				case 3:
					return getEquatorialFaceCoords(0, i, j);
					break;
				case 4:
					return getEquatorialFaceCoords(90, i, j);
					break;
				case 5:
					return getPolarFaceCoords(false, i, j);
					break;
				default:
					std::cerr << "\033[1;31mError! Incorrect face index \033[1;37m" << face << "\033[1;31m!\033[0m" << std::endl;
					return Coordinate(0, 0);
			}
		}

		template<typename T>
		const Coordinate cubedSphereDS<T>::getPolarFaceCoords(const bool NP, const size_t i, const size_t j) const{
			#if __cplusplus > 201703L
				double pi = std::numbers::pi;
			#else
				double pi = 3.1415926535;
			#endif
			double edgeLength = PWM::Utils::radToDeg((pi/2) / this->getFaceWidth());
			double x = ((double) j + 0.5) - this->getFaceWidth() / 2;
			double y = -(((double) i + 0.5) - this->getFaceWidth() / 2);
			double longi = 45 + (NP ? PWM::Utils::radToDeg(std::atan2(y, x)) : PWM::Utils::radToDeg(-std::atan2(y, x)));
			
			if (longi > 180)
				longi -= 360;
			else if (longi < -180)
				longi += 360;
			else if (longi == 180)
				longi = -180;

			double lat = 90 - ((std::max(std::abs(x), std::abs(y))) * edgeLength);
			if (!NP)
				lat = -lat;
			
			return Coordinate(lat, longi);
		}

		template<typename T>
		const Coordinate cubedSphereDS<T>::getEquatorialFaceCoords(const double leftEdgeLongitude, const size_t i, const size_t j) const{
			double lat = 0;
			double longi = 0;
			int fw = this->getFaceWidth();
			double edgeLength = 90. / fw;
			if (fw % 2){
				int latIndex = std::abs((fw / 2.) - ((int) i));
				lat = latIndex * edgeLength;
				lat = (i < (fw / 2)) ? lat : -lat;

				longi = leftEdgeLongitude + j * edgeLength;
			}else{
				lat = -((edgeLength * (int) i) + (edgeLength / 2) - 45.);
				longi = leftEdgeLongitude + j * edgeLength + edgeLength / 2.;
			}
			return Coordinate(lat, longi);
		}
		
		template<typename T>
		const size_t cubedSphereDS<T>::convertFaceCoordToIndex(const size_t face, const size_t i, const size_t j) const{
			return (face * getFaceWidth() * getFaceWidth()) + (i * getFaceWidth()) + j;
		}

		template<typename T>
		const std::tuple<size_t, size_t, size_t> cubedSphereDS<T>::convertIndexToFaceCoord(const size_t index) const{
			size_t face, i, j;
			face = index / (getFaceWidth() * getFaceWidth());
			i = (index - (face * getFaceWidth() * getFaceWidth())) / getFaceWidth();
			j = index - (face * getFaceWidth() * getFaceWidth()) - (i * getFaceWidth());
			return std::make_tuple(face, i, j);
		}

		template<typename T>
		void swap(cubedSphereDS<T>& i, cubedSphereDS<T>& j){
			std::swap(i.faceWidth, j.faceWidth);
			std::swap((AbstractPWMDataStructure<T>) i, (AbstractPWMDataStructure<T>) j);
		}
    }
}

#endif //CUBED_SPHERE_DATA_STRUCTURE_H