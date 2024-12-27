/*******************************************************************************
 *
 * EcoSynth - Data-driven Authoring of Large-Scale Ecosystems (Undergrowth simulator)
 * Copyright (C) 2022  J.E. Gain  (jgain@cs.uct.ac.za)
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


#include "noisefield.h"

NoiseField::NoiseField(Terrain * ter, int dstep, long sval)
{
    int tx, ty;

    terrain = ter;
    terrain->getGridDim(tx, ty);
    dimx = tx * dstep;
    dimy = tx * dstep;
    nmap = new MapFloat();
    nmap->setDim(dimx, dimy);
    nmap->fill(0.0f);
    seed = sval;
    dice = new DiceRoller(0, 1000, sval);
    init();
}

void NoiseField::init()
{
    for(int x = 0; x < dimx; x++)
        for(int y = 0; y < dimy; y++)
        {
            nmap->set(x, y, (float) dice->generate() / 1000.0f);
        }
}

float NoiseField::getNoise(vpPoint p)
{
    float tx, ty, convx, convy;
    int x, y;

    terrain->getTerrainDim(tx, ty);
    convx = (float) (dimx-1) / tx;
    convy = (float) (dimy-1) / ty;
    x = (int) (p.x * convx);
    y = (int) (p.z * convy);
    return nmap->get(x, y);
}
