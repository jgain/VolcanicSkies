/*******************************************************************************
 *
 * EcoSynth - Data-driven Authoring of Large-Scale Ecosystems (Undergrowth simulator)
 * Copyright (C) 2020  J.E. Gain  (jgain@cs.uct.ac.za)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 ********************************************************************************/


#ifndef NOISEFIELD_H
#define NOISEFIELD_H

#include "terrain.h"
#include "dice_roller.h"
#include "common/basic_types.h"

class NoiseField
{
private:
    long seed;      //< random number initialization to ensure a Noisefield can be reproduced
    Terrain * terrain;  //< underlying terrain
    int dimx, dimy; //< dimsions of the noisefield (a multiplicative factor of the terrain)
    MapFloat * nmap;  //< field of noise values
    DiceRoller * dice; //< random number generator

public:

    /**
     * @brief NoiseField initializer
     * @param ter   underlying terrain
     * @param dstep number of sub-cells per terrain cell (basically a multiplication factor)
     * @param sval  seed value
     */
    NoiseField(Terrain * ter, int dstep, long sval);

    ~NoiseField(){ delete dice; }

    /// initialize the noisemap
    void init();

    /// recover a random value in [0, 1] at a point on the terrain
    float getNoise(vpPoint p);
};

#endif // NOISEFIELD_H
