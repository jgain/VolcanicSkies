#ifndef PWM_MODEL_SUN_H
#define PWM_MODEL_SUN_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include "json.hpp"
namespace PWM{
    namespace Model{
        template<typename T> class sun;

        template<typename T>
        void to_json(nlohmann::json& j, const sun<T>& s);
        
        template<typename T>
        void from_json(nlohmann::json& j, const sun<T>& s);
        
        template<typename T> class sun{
            private:
                /**
                 * The name of the star/sun
                 */
                std::string name;
                
                /**
                 * How much power the sun/star puts out ignoring the atmospheric or other absorbtion
                 * in Watt per meter squared at a distance of 1 AU.
                 */
                T powerOutput;
                
                /**
                 * The sun/star's position in x,y,z coordinates.
                 *
                PWM::PWMDataStructure::Position pos;*/

                /**
                 * The distance from the sun to the planet in AU.
                 */
                T distance;                
                /**
                 * The wavelength at which this sun/star sends out the majority of its energy
                 * i.e., the peak of the black body radiation graph. Measured in nanometres.
                 * Not used for now, can be an expansion if the atmosphere absorption/opacity
                 * is put in. 550 for the Sun.
                 */
                T colour;

                /**
                 * The apparent declination of the sun, used for ground heating
                 * determination w.r.t. the seasons. Varies between ± 23.5 deg
                 * for Earth.
                 */
                T apparentDeclination;

                /**
                 * The direction in which the seasons are changing. +1 if the northern hemisphere
                 * of the planet is moving to summer, -1 if the southern hemisphere is
                 * moving to summer.
                 */
                int seasonDirection;

                /**
                 * The apparent right ascension of this sun, used for ground heating
                 * determination w.r.t. the day night cycle. Varies between 0 and 360
                 * deg for Earth.
                 */
                T apparentRightAscension;

                /**
                 * The planet's axial tilt with reference to its orbit around this
                 * sun. 23.5 deg for Earth.
                 */
                T axialTilt;

                /**
                 * The period of the planet's axial tilt w.r.t. this sun.
                 * On Earth, equivalent to solar year (not sidereal year),
                 * i.e, 31 556 923.488 seconds.
                 */
                T tropicalYearLength;

                /**
                 * The solar day period for the planet with reference to this sun
                 * in seconds. 86400 for Earth.
                 */
                T solarDayLength;
            public:
                sun();
                sun(std::string JSONfile);
                int writeSunFile(std::string file) const;

                const std::string getName() const;
                const T& getPower() const;
                const T& getColour() const;
                //const T& getPosition() const;
                const T& getDistance() const;
                const T& getApparentDeclination() const;
                const int getSeasonDirection() const;
                const T& getApparentRightAscension() const;
                const T& getAxialTilt() const;
                const T& getTropicalYearLength() const;
                const T& getSolarDayLength() const;

                void setName(std::string n);
                void setApparentDeclination(T dec);
                void setSeasonDirection(int SD);
                void flipSeasonDirection();
                void setApparentRightAscension(T ra);
                void setPower(T pow);                
                void setColour(T col);
                void setDistance(T dis);
                void setAxialTilt(T axt);
                void setTropicalYearLength(T tyl);
                void setSolarDayLength(T sdl);

                bool operator==(const sun& rhs) const;
                bool operator!=(const sun& rhs) const;

                template<typename U> operator sun<U>(){
                    sun<U> res = sun<U>();
                    res.setPower((U) this->powerOutput);
                    res.setColour((U) this->colour);
                    res.setDistance((U) this->distance);
                    res.setApparentDeclination((U) this->apparentDeclination);
                    res.setSeasonDirection(this->seasonDirection);
                    res.setApparentRightAscension((U) this->apparentRightAscension);
                    res.setAxialTilt((U) this->axialTilt);
                    res.setTropicalYearLength((U) this->tropicalYearLength);
                    res.setSolarDayLength((U) this->solarDayLength);

                    return res;
                }
        };

        template<typename T>
    inline sun<T>::sun(){
            setApparentDeclination(0);
            setApparentRightAscension(0);
            seasonDirection = 1;
        }
        
        template<typename T>
    inline sun<T>::sun(std::string JSONfile){
            std::ifstream fil(JSONfile);
            if (!fil.good()){
                std::cout << "\033[1;31mError! Planet json file " << JSONfile << " not found!\033[0m" << std::endl;
                return;
            }
            nlohmann::json jsonData;
            try{
                fil >> jsonData;
                *this = jsonData.get<sun>();
            }
            catch (nlohmann::json::exception& e){
                std::cout << "\033[1;31mError! Sun json file not in correct format!\033[0m" << std::endl;
                std::cout << "\033[1;37mSun json file not read in.\033[0m" << std::endl;
            }
            setApparentDeclination(0);
            setApparentRightAscension(0);
            seasonDirection = 1;
        }

        template<typename T>
    inline int sun<T>::writeSunFile(std::string file) const{
            std::ofstream fil(file);
            if (!fil.good()){
                std::cerr << "\033[1;31mError! Sun file " << file << " unsuitable for writing!\033[0m" << std::endl;
                return -1;
            }
            try {
                nlohmann::json jsonData = *this;
                fil << std::setw(4) << jsonData << std::endl;
                return 0;
            } catch (nlohmann::json::exception& e) {
                std::cerr << "\033[1;31mError when writing sun json file!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mSun json file not written.\033[0m" << std::endl;
                return -2;
            }
        }

        template<typename T>
    inline const std::string sun<T>::getName() const{
            return name;
        }

        template<typename T>
    inline const T& sun<T>::getPower() const{
            return powerOutput;
        }

        template<typename T>
    inline const T& sun<T>::getColour() const{
            return colour;
        }

        template<typename T>
    inline const T& sun<T>::getDistance() const{
            return distance;
        }

        template<typename T>
    inline const T& sun<T>::getApparentDeclination() const{
            return apparentDeclination;
        }

        template<typename T>
    inline const int sun<T>::getSeasonDirection() const{
            return seasonDirection;
        }

        template<typename T>
    inline const T& sun<T>::getApparentRightAscension() const{
            return apparentRightAscension;
        }

        template<typename T>
    inline const T& sun<T>::getTropicalYearLength() const{
            return tropicalYearLength;
        }

        template<typename T>
    inline const T& sun<T>::getAxialTilt() const{
            return axialTilt;
        }
        
        template<typename T>
    inline const T& sun<T>::getSolarDayLength() const{
            return solarDayLength;
        }

        template<typename T>
    inline void sun<T>::setName(std::string n){
            name = n;
        }
        
        template<typename T>
    inline void sun<T>::setApparentDeclination(T dec){
            apparentDeclination = dec;
        }

        template<typename T>
    inline void sun<T>::setSeasonDirection(int SD){
            seasonDirection = SD;
        }

        template<typename T>
    inline void sun<T>::flipSeasonDirection(){
            seasonDirection *= -1;
        }

        template<typename T>
    inline void sun<T>::setApparentRightAscension(T ra){
            apparentRightAscension = ra;
        }

        template<typename T>
    inline void sun<T>::setPower(T pow){
            powerOutput = pow;
        }

        template<typename T>
    inline void sun<T>::setColour(T col){
            colour = col;
        }

        template<typename T>
    inline void sun<T>::setDistance(T dis){
            distance = dis;
        }

        template<typename T>
    inline void sun<T>::setAxialTilt(T axt){
            axialTilt = axt;
        }

        template<typename T>
    inline void sun<T>::setTropicalYearLength(T tyl){
            tropicalYearLength = tyl;
        }

        template<typename T>
    inline void sun<T>::setSolarDayLength(T sdl){
            solarDayLength = sdl;
        }

        template<typename T>
    inline bool sun<T>::operator==(const sun<T>& rhs) const{
            if (this->getPower() != rhs.getPower())
                return false;
            if (this->getColour() != rhs.getColour())
                return false;
            if (this->getDistance() != rhs.getDistance())
                return false;
            if (this->getApparentDeclination() != rhs.getApparentDeclination())
                return false;
            if (this->getSeasonDirection() != rhs.getSeasonDirection())
                return false;
            if (this->getApparentRightAscension() != rhs.getApparentRightAscension())
                return false;
            if (this->getAxialTilt() != rhs.getAxialTilt())
                return false;
            if (this->getTropicalYearLength() != rhs.getTropicalYearLength())
                return false;
            if (this->getSolarDayLength() != rhs.getSolarDayLength())
                return false;
            return true;
        }

        template<typename T>
    inline bool sun<T>::operator!=(const sun<T>& rhs) const{
            return !(*this == rhs);
        }

        //JSON lib interaction functions
        template<typename T>
    inline void to_json(nlohmann::json& j, const sun<T>& s){
            j = nlohmann::json{
                {"name", s.getName()},
                {"powerOutput", s.getPower()},
                {"colour", s.getColour()},
                {"distance", s.getDistance()},
                {"axialTilt", s.getAxialTilt()},
                {"tropicalYearLength", s.getTropicalYearLength()},
                {"solarDayLength", s.getSolarDayLength()}
            };
        }

        template<typename T>
    inline void from_json(const nlohmann::json& j, sun<T>& s){
            s.setName(j.at("name").get<std::string>());
            s.setPower(j.at("powerOutput").get<T>());
            s.setColour(j.at("colour").get<T>());
            s.setDistance(j.at("distance").get<T>());
            s.setAxialTilt(j.at("axialTilt").get<T>());
            s.setTropicalYearLength(j.at("tropicalYearLength").get<T>());
            s.setSolarDayLength(j.at("solarDayLength").get<T>());
        }
    }
}

#endif //PWM_MODEL_SUN_H
