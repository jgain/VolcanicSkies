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
#include <chrono>
#include "flatStaggeredGrid.h"
#include <iostream>
#include <memory>
#include <random>

typedef double valType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;
using namespace PWM::PWMDataStructure;

int main(int argc, char *argv[])
{
    auto source = std::make_shared<dsType>(dsType(1024, 1024));
    auto dest = std::make_shared<dsType>(dsType(1024, 1024));
    auto dest2 = std::make_shared<dsType>(dsType(1024, 1024));
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(-100,100); // distribution in range [1, 6]
    for (int i = 0; i < source->size(); ++i){
        source->setData(i, dist(rng));
    }
    source->movingAverageSmoothingA(dest, 13);
    source->movingAverageSmoothingWFA(dest2, 13);
    if (*dest != *dest2){
        std::cerr << "\033[1;31mWelford's Algorithm smoothing is not equal!\033[0m" << std::endl;
        return -1;
    }

    double timeSpent = 0;
    int loops = 10;
    std::chrono::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::time_point t2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < loops; ++i){
        dest->copy(0);
        t1 = std::chrono::high_resolution_clock::now();
        source->movingAverageSmoothingA(dest, 12);
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> msDouble = t2 - t1;
        timeSpent += msDouble.count();
    }
    std::cout << "Average time over 10 runs for smoothing with outermost loop parallelised: " << timeSpent / loops << " milliseconds." << std::endl;

    timeSpent = 0;
    for (int i = 0; i < loops; ++i){
        dest->copy(0);
        t1 = std::chrono::high_resolution_clock::now();
        source->movingAverageSmoothingB(dest, 12);
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> msDouble = t2 - t1;
        timeSpent += msDouble.count();
    }
    std::cout << "Average time over 10 runs for smoothing with inner loop parallelised: " << timeSpent / loops << " milliseconds." << std::endl;

    timeSpent = 0;
    for (int i = 0; i < loops; ++i){
        dest->copy(0);
        t1 = std::chrono::high_resolution_clock::now();
        source->movingAverageSmoothingC(dest, 12);
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> msDouble = t2 - t1;
        timeSpent += msDouble.count();
    }
    std::cout << "Average time over 10 runs for smoothing with summation loop parallelised: " << timeSpent / loops << " milliseconds." << std::endl;

    timeSpent = 0;
    for (int i = 0; i < loops; ++i){
        dest->copy(0);
        t1 = std::chrono::high_resolution_clock::now();
        source->movingAverageSmoothingWFA(dest, 12);
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> msDouble = t2 - t1;
        timeSpent += msDouble.count();
    }
    std::cout << "Average time over 10 runs for smoothing with parallel Welford algorithm: " << timeSpent / loops << " milliseconds." << std::endl;
    return 0;
}

