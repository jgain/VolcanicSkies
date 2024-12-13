#ifndef PWM_SSG_DS
#define PWM_SSG_DS
// This class serves as a class template for the staggered spherical grid data structure.
// Cilliers Pretorius
// 05 May 2021
#include <cmath>
#include <cstddef>
#include <iostream>
#if __cplusplus > 201703L
	#include <numbers>
#endif
#include <omp.h>
#include <random>
#include <tuple>
#include <typeinfo>
#include <utility>
#include <ctime>
#include "AbstractPWMDataStructure.h"
#include "mathUtils.h"

namespace PWM {
	namespace PWMDataStructure{
		template<typename T>
		class SSGDataStructure;
		
		template<typename T>
		void swap(SSGDataStructure<T>& i, SSGDataStructure<T>& j);

		template<typename T>
		class SSGDataStructure : public AbstractPWMDataStructure<T> {
			private:
				size_t theta; //latitude
				size_t phi; //longitude
			/**
			 * Declarations for all the functions of this class template.
			 */
			public:
				//Default constructor
				inline SSGDataStructure<T>(size_t theta, size_t phi);//done
				
				/*//Default destructor
				inline ~SSGDataStructure();//done
			
				//Default copy constructor
				inline SSGDataStructure(const SSGDataStructure<T> & other);//done
			
				//Default move constructor
				inline SSGDataStructure(SSGDataStructure<T> && other);//done
			
				//Default copy assignment
				inline SSGDataStructure<T> & operator= (const SSGDataStructure<T> & other);//done
			
				//Default move assignment
				inline SSGDataStructure<T> & operator= (SSGDataStructure<T> && other);//done*/
				
				//Equality comparator
				inline bool operator== (const SSGDataStructure<T> & other) const;//done
				inline bool operator!= (const SSGDataStructure<T> & other) const;//done
			
				//safe read
				inline const T getData(const size_t index) const override;//done
				
				//Read access for i, j calls, meant for the advection engine.
				inline const T getData(const int i, const int j) const;
				
				//Utility function that bilinearly interpolates values to find the center
				inline const T getLerpedVal(const int index) const;
				inline const T getLerpedVal(const double i, const double j) const; 

				//Function that randomly initialises the datastructure with the given min and max values
				void randomInit(const T& min, const T& max);
				//Get the sum of all elements in the data structure.
				inline const T sum() const;
			
				//Get the average of all elements in the data structure.
				inline const T mean() const;
			
				//safe write to data point
				inline void setData(const int index, const T& val);//done
				
				//Write access for i, j calls, meant for the advection and convection engines.
				inline void setData(const int i, const int j, const T& val);//done
				
				//Retrieve coordinates for a specific point (used wherever the spherical surface influences equation)
				inline Coordinate getCoordinates(const size_t index) const override;//done
				inline Coordinate getCoordinates(const size_t i, const size_t j) const;//done
			
				inline size_t getPhi() const;//done
				inline size_t getTheta() const;//done
				inline const size_t size() const;//done
				inline const double gridLength() const;
				inline void copy(const T& val);//done
				inline SSGDataStructure<T> & operator= (const T& val);//done
				inline void copy(const SSGDataStructure<T> & other);//done
				
				//small utility function to check if two structures are of equal size
				inline bool checkSize(const SSGDataStructure<T> & other) const;//done
				
				//small utility function to convert i,j reference to linear access index
				inline int convert2Dto1D(size_t i, size_t j) const;//done
				inline int convert2Dto1D(int i, int j) const;//done

				//small utility function to convert index reference to i, j references
				inline std::pair<int, int> convert1Dto2D(size_t index) const;
				inline std::pair<int, int> convert1Dto2D(int index) const;

				friend void swap <> (SSGDataStructure<T>& i, SSGDataStructure<T>& j);
		};
		
		/**
		 * Implementations of all the functions declared above.
		 */
		template<typename T>
		SSGDataStructure<T>::SSGDataStructure(size_t theta, size_t phi): AbstractPWMDataStructure<T>(theta * phi){
			this->theta = theta;
			this->phi = phi;
		}
		
		template<typename T>
		const size_t SSGDataStructure<T>::size() const{
			return getTheta() * getPhi();
		}

		/*template<typename T>
		SSGDataStructure<T>::~SSGDataStructure<T>(){
		}

		template<typename T>
		SSGDataStructure<T>::SSGDataStructure(const SSGDataStructure<T> & other): AbstractPWMDataStructure<T>(theta * phi), theta(other.getTheta()), phi(other.getPhi()){
			copy(other);
		}
		
		template<typename T>
		SSGDataStructure<T>::SSGDataStructure(SSGDataStructure<T> && other): AbstractPWMDataStructure<T>(theta * phi), theta(other.getTheta()), phi(other.getPhi()){
			this->data = std::move(other.data);
		}
		
		template<typename T>
		SSGDataStructure<T> & SSGDataStructure<T>::operator= (const SSGDataStructure<T> & other) {
			if (this == other)
				return *this;
			theta = other.getTheta();
			phi = other.getPhi();
			AbstractPWMDataStructure<T>(theta * phi);
			copy(other);
			
			return *this;
		}
		
		template<typename T>
		SSGDataStructure<T> & SSGDataStructure<T>::operator= (SSGDataStructure<T> && other) {
			if (this == other)
				return *this;
			theta = other.getTheta();
			phi = other.getPhi();
			AbstractPWMDataStructure<T>(theta * phi);
			copy(other);
			
			other->data = nullptr;
			return *this;
		}*/
		
		template<typename T>
		bool SSGDataStructure<T>::operator== (const SSGDataStructure<T> & other) const {
			if (typeid(getData(0)) != typeid(other.getData(0))){
				return false;
			}
			if (!checkSize(other)){
				return false;
			}
			for (int i = 0; i < size(); ++i){
				if (getData(i) != other.getData(i)) {
					return false;
				}
			}
			return true;
		}
		
		template<typename T>
		bool SSGDataStructure<T>::operator!= (const SSGDataStructure<T> & other) const {
			return !(*this == other);
		}
		
		template<typename T>
		bool SSGDataStructure<T>::checkSize(const SSGDataStructure<T> & other) const {
			return (getTheta() == other.getTheta()) && (getPhi() == other.getPhi());
		}
		
		template<typename T>
		const T SSGDataStructure<T>::getData(const size_t index) const {
			if (index < size())
				return this->data[index];
			std::cerr << "Error: Out of bounds array read access in SSGDataStructure::getData(int index) with index = " << index << "!" << std::endl;
			return T();
		}
		
		template<typename T>
		const T SSGDataStructure<T>::getData(const int i, const int j) const {
			return getData(convert2Dto1D(i, j));
		}
		
		template<typename T>
		const T SSGDataStructure<T>::getLerpedVal(const double i, const double j) const {
			//Get theta/phi indices
			auto normalizedTheta = i * (1.0 / gridLength());
			auto normalizedPhi = j * (1.0 / gridLength());

			int thetaIndex = static_cast<int>(std::floor(normalizedTheta));
			int phiIndex = static_cast<int>(std::floor(normalizedPhi));

			auto alphaTheta = normalizedTheta - static_cast<T>(thetaIndex);
			auto alphaPhi = normalizedPhi - static_cast<T>(phiIndex);

			size_t thetaLower = thetaIndex;
			size_t thetaUpper = thetaIndex + 1;
			size_t phiLower = phiIndex % getPhi();
			size_t phiUpper = (phiIndex + 1) % getPhi();

			auto lower = PWM::Utils::interpolate<T>(getData(thetaLower, phiLower), getData(thetaLower, phiUpper), alphaPhi);
			auto upper = PWM::Utils::interpolate<T>(getData(thetaUpper, phiLower), getData(thetaUpper, phiUpper), alphaPhi);
			
			return PWM::Utils::interpolate<T>(lower, upper, alphaTheta);
		}

		template<typename T>
		const T SSGDataStructure<T>::getLerpedVal(const int index) const {
			auto c = convert1Dto2D(index);
			return getLerpedVal(c.first, c.second);
		}

		template<typename T>
		void SSGDataStructure<T>::setData(const int index, const T& val) {
			if (index < size()){
				this->data[index] = val;
				return;
			}
			std::cerr << "Error: Out of bounds array write access in SSGDataStructure::getData(int index) with index = " << index << "!" << std::endl;
			return;
		}
		
		template<typename T>
		void SSGDataStructure<T>::setData(const int i, const int j, const T& val) {
			setData(convert2Dto1D(i, j), val);
			return;
		}
		
		template<typename T>
		size_t SSGDataStructure<T>::getTheta() const {
			return this->theta;
		}
		
		template<typename T>
		size_t SSGDataStructure<T>::getPhi() const {
			return this->phi;
		}
		
		template<typename T>
		int SSGDataStructure<T>::convert2Dto1D(size_t i, size_t j) const {
			if (i < 0)
				i = getTheta() - std::abs((int) i);
			if (j < 0)
				j = getPhi() - std::abs((int) j);
			return ((i % getTheta())  * getPhi()) + (j % getPhi());
		}

		template<typename T>
		int SSGDataStructure<T>::convert2Dto1D(int i, int j) const{
			if (i < 0)
				i = getTheta() - std::abs(i);
			if (j < 0)
				j = getPhi() - std::abs(j);
			return ((i % getTheta())  * getPhi()) + (j % getPhi());
		}

		template<typename T>
		std::pair<int, int> SSGDataStructure<T>::convert1Dto2D(size_t index) const {
			int i, j;
			if (index < 0)
				index = (size() - std::abs((int) index));
			i = index / getPhi();
			j = index % getPhi();
			return std::make_pair(i, j);
		}

		template<typename T>
		std::pair<int, int> SSGDataStructure<T>::convert1Dto2D(int index) const{
			int i, j;
			if (index < 0)
				index = size() - std::abs(index);
			i = index / getPhi();
			j = index % getPhi();
			return std::make_pair(i, j);
		}
		
		template<typename T>
		void SSGDataStructure<T>::copy(const T& val) {
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i)
					this->data[i] = val;
		}
		
		template<typename T>
		SSGDataStructure<T> & SSGDataStructure<T>::operator= (const T& val) {
			copy(val);
			return *this;
		}
		
		template<typename T>
		void SSGDataStructure<T>::copy(const SSGDataStructure<T> & other) {
			if (!checkSize(other))
				std::cerr << "Error: Data structures of unequal sizes in SSGDataStructure::copy(const SSGDataStructure<T> & other)!" << std::endl;
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i)
					this->data[i] = other.getData(i);
            #pragma omp barrier
		}
		
		template<typename T>
		Coordinate SSGDataStructure<T>::getCoordinates(const size_t index) const{
			auto r = convert1Dto2D(index);
			return getCoordinates(r.first, r.second);
		}

		template<typename T>
		Coordinate SSGDataStructure<T>::getCoordinates(const size_t i, const size_t j) const{
			double lat = 0;
			double longi = 0;
			if (getTheta() % 2){
				size_t latIndex = std::abs((int) ((getTheta() / 2) - i));
				lat = latIndex * (90 / (getTheta() / 2));
				lat = (i < (getTheta() / 2)) ? lat : -lat;
			}else{
				lat = -(((180 / (double) getTheta()) * (double) i) + (90 / (double) getTheta()) - 90);//Theta-0 is north pole, i.e. 90 deg.
			}
			if (getPhi() % 2){
				size_t longiIndex = std::abs((int) ((getPhi() / 2) - j));
				longi = longiIndex * (180 / (getPhi() / 2));
				longi = (j < (getPhi() / 2)) ? -longi : longi;
			}else{
				longi = ((360 / (double) getPhi()) * (double) j) + (180 / (double) getPhi()) - 180;//Phi-0 is 180 deg west, i.e., -180 deg.
			}
			auto res = Coordinate(lat, longi);
			return res;
		}
		
		template<typename T>
		const T SSGDataStructure<T>::sum() const{
			auto res = T();
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i)
					#pragma omp atomic
						res += this->data[i];
			#pragma omp barrier
			return res;
		}
		
		template<typename T>
		const T SSGDataStructure<T>::mean() const{
			return sum() / size();
		}

		template<typename T>
		const double SSGDataStructure<T>::gridLength() const {
			#if __cplusplus > 201703L
				return std::numbers::pi / getTheta();
			#else
				return 3.1415926535 / getTheta();
			#endif
		}

		template<typename T>
		void SSGDataStructure<T>::randomInit(const T& min, const T& max){
			std::mt19937_64 mt(time(NULL));
			auto d = safeMinMax(min, max);
			std::uniform_real_distribution<T> dist(std::get<0>(d), std::get<1>(d));
			for (int i = 0; i < size(); ++i)
				setData(i, (dist(mt) - std::get<2>(d)));
		}

		template<typename T>
		std::tuple<T, T, T> safeMinMax(T& min, T& max){
			if (min < 0){
				return std::make_tuple((T) 0, max + std::abs(min), std::abs(min));
			}
			return std::make_tuple(min, max, (T) 0);
		}

		template<typename T>
		void swap(SSGDataStructure<T>& i, SSGDataStructure<T>& j){
			std::swap(i.data, j.data);
		}
	}
}

#endif //PWM_SSG_DS
