#ifndef PWM_PWMDATASTRUCTURE_POSITION_H
#define PWM_PWMDATASTRUCTURE_POSITION_H

#include <tuple>
namespace PWM{
    namespace PWMDataStructure{
        class Position{
            private:
                double x, y, z;
            public:
                Position(double X = 0, double Y = 0, double Z = 0);
                
                const double getX() const;
                const double getY() const;
                const double getZ() const;
                
                const double getDistanceTo(Position& other) const;
                
                void move(double newX, double newY, double newZ);
                void move(std::tuple<double, double, double>& distance);
                void move(std::tuple<double, double, double>& velocity, double dt);
                void move(Position& newPos);

                bool operator==(Position& other);
                bool operator!=(Position& other);

        };
    }
}
#endif //PWM_PWMDATASTRUCTURE_POSITION_H