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


// volcano.cpp: for managing volcano and weather data
// author: James Gain
// date: June 2023

#ifdef _WIN32
#include <glew.h>
#else
#include <GL/glew.h>
#endif

#include "volcano.h"
#include <sstream>
#include <streambuf>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <sys/stat.h>
#include <QImage>

using namespace std;

QImage & AirLayers::getCurrentLayerImage()
{
    switch(vmode)
    {
        case VolcanoMode::TEMPERATURE:
            return getCurrentTemperature();
            break;
        case VolcanoMode::MOISTURE:
            return getCurrentMoisture();
            break;
        case VolcanoMode::CONVECTION:
            return getCurrentConvection();
            break;
        case VolcanoMode::ASH:
            return getCurrentAsh();
            break;
        case VolcanoMode::VELOCITYX:
            return getCurrentVelocityX();
            break;
        case VolcanoMode::VELOCITYY:
            return getCurrentVelocityY();
            break;
        case VolcanoMode::CONDENSED:
            return getCurrentCondensed();
            break;
        //case VolcanoMode::VMEND:
        //    return nullptr;
        //    break;
    }
}

void AirLayers::loadAirLayers(string &basefile, int interval)
{
    ifstream infile;
    string filename;
    
    clear();
    
    filename = basefile + "layerstruct.txt";
    infile.open((char *) filename.c_str(), ios_base::in);
    if(infile.is_open())
    {
        // load airlayers structure
        infile >> numlayers;
        
        for(int a = 0; a < numlayers; a++)
        {
            AirLayer plane;
            infile >> plane.altitude;
            // cerr << plane.altitude << endl;
            plane.active = false;
            planes.push_back(plane);
        }
        infile.close();
    }
    else
    {
        cerr << "Error AirLayers::loadAirLayers - unable to open file " << filename << endl;
    }
    if(numlayers > 0)
        setActiveStatus(true);
    
    /*
    // load temperature
    for(int f = 0; f < numframes; f++)
    {
        std::vector<QImage> t;
        
        for(int a = 0; a < numlayers; a++)
        {
            QImage inimg, outimg;
            int alt = planes[a].altitude;
            filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_Temperature_Step_" + to_string((f+1)*50) + ".png";
           
            // load image
            if(!inimg.load(QString::fromStdString(filename)))
                cerr << "Error AirLayers::loadAirLayers - unable to open file " << filename << endl;
            
            // Qt prep image for OpenGL
            inimg = inimg.convertToFormat(QImage::Format_ARGB32);
            outimg = QGLWidget::convertToGLFormat(inimg);
            t.push_back(outimg);
        }
        temperature.push_back(t);
    }
    
    // load moisture
    for(int f = 0; f < numframes; f++)
    {
        std::vector<QImage> m;
        
        for(int a = 0; a < numlayers; a++)
        {
            QImage inimg, outimg;
            int alt = planes[a].altitude;
            filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_Moisture_Step_" + to_string((f+1)*interval) + ".png";
           
            // load image
            if(!inimg.load(QString::fromStdString(filename)))
                cerr << "Error AirLayers::loadAirLayers - unable to open file " << filename << endl;

            // Qt prep image for OpenGL
            QImage fixedimg(inimg.width(), inimg.height(), QImage::Format_RGB32);
            outimg = QGLWidget::convertToGLFormat(fixedimg);
            m.push_back(outimg);
        }
        moisture.push_back(m);
    }
    
    // load convection
    for(int f = 0; f < numframes; f++)
    {
        std::vector<QImage> c;
        
        for(int a = 0; a < numlayers-1; a++)
        {
            QImage inimg, outimg;
            filename = basefile + "ConvectionLayer_" + to_string(a) + "_Vertical_Velocity_Step_" + to_string((f+1)*interval) + ".png";
           
            // load image
            if(!inimg.load(QString::fromStdString(filename)))
                cerr << "Error AirLayers::loadAirLayers - unable to open file " << filename << endl;

            // Qt prep image for OpenGL
            QImage fixedimg(inimg.width(), inimg.height(), QImage::Format_RGB32);
            outimg = QGLWidget::convertToGLFormat(fixedimg);
            c.push_back(outimg);
        }
        convection.push_back(c);
    }*/
    
    // load images
    for(auto v: all_volcanomodes)
        for(int f = 0; f < numframes; f++)
        {
            std::vector<QImage> imgs;
            int currlayers = numlayers;
            if(v == VolcanoMode::CONVECTION)
                currlayers -= 1;
            
            for(int lyr = 0; lyr < currlayers; lyr++)
            {
                int a = lyr*2;
                QImage inimg, outimg;
                int alt = planes[lyr].altitude;
                filename = "";
                
                switch(v)
                {
                    case VolcanoMode::TEMPERATURE:
                        filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_Temperature_Step_" + to_string(f*interval) + ".png";
                        break;
                    case VolcanoMode::MOISTURE:
                        /*
                        filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_Moisture_Step_" + to_string(f*interval) + ".png";*/
                        break;
                    case VolcanoMode::CONVECTION:
                        filename = basefile + "ConvectionLayer_" + to_string(a) + "_Vertical_Velocity_Step_" + to_string(f*interval) + ".png";
                        break;
                    case VolcanoMode::ASH:
                        filename = basefile + "Layer_" + to_string(a) + "_(" +
                        to_string(alt) + "m)_Particulates_Step_" + to_string(f*interval) + ".png";
                        break;
                    case VolcanoMode::VELOCITYX:
                        /*
                        filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_VelocityX_Step_" + to_string(f*interval) + ".png";*/
                        break;
                    case VolcanoMode::VELOCITYY:
                        /*
                        filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_VelocityY_Step_" + to_string(f*interval) + ".png";*/
                        break;
                    case VolcanoMode::CONDENSED:
                        filename = basefile + "Layer_" + to_string(a) + "_(" + to_string(alt) + "m)_CondensedWater_Step_" + to_string(f*interval) + ".png";
                        break;
                    case VolcanoMode::VMEND:
                        break;
                }
              

                if(filename != "")
                {
                    // load image
                    if(!inimg.load(QString::fromStdString(filename), "png"))
                        cerr << "Error AirLayers::loadAirLayers - unable to open file " << filename << endl;

                    // Qt prep image for OpenGL
                    inimg = inimg.convertToFormat(QImage::Format_ARGB32);
                    outimg = QGLWidget::convertToGLFormat(inimg);
                    
                    imgs.push_back(outimg);
                }
            }
            switch(v)
            {
                case VolcanoMode::TEMPERATURE:
                    temperature.push_back(imgs);
                    break;
                case VolcanoMode::MOISTURE:
                    moisture.push_back(imgs);
                    break;
                case VolcanoMode::CONVECTION:
                    convection.push_back(imgs);
                    break;
                case VolcanoMode::ASH:
                    ash.push_back(imgs);
                    break;
                case VolcanoMode::VELOCITYX:
                    velocityx.push_back(imgs);
                    break;
                case VolcanoMode::VELOCITYY:
                    velocityy.push_back(imgs);
                    break;
                case VolcanoMode::CONDENSED:
                    condensed.push_back(imgs);
                    break;
                case VolcanoMode::VMEND:
                    break;
            }
           
        }
}

/*
void Skirts::getGridDim(int & dx, int & dy) const
{
    // assumes all grids have the same size
    if((int) skirtlayers.size() > 0)
    {
        if((int) skirtlayers[0].heights.size() > 0)
        {
            dx = skirtlayers[0].heights[0]->width();
            dy = skirtlayers[0].heights[0]->height();
        }
        else
        {
            cerr << "Error skirt::getGridDim - trying to access an empty grid" << endl;
        }
    }
}*/

int Skirts::getNumSamples(int l)
{
    return (int) skirtlayers[l].frames[currframe].size();
}

void Skirts::getFrameSample(int l, int s, vpPoint &pos, float &cval)
{
    pos = skirtlayers[l].frames[currframe][s].position;
    cval = skirtlayers[l].frames[currframe][s].condensation;
    
    /*
    float hght = skirtlayers[l].heights[currframe]->get(x, y);
    cval = skirtlayers[l].condensations[currframe]->get(x, y);
    
    // find position in world coordinate space
    int gx, gy;
    float convx, convy;

    getGridDim(gx, gy);
    convx = skirtlayers[l].extent.i / (float) (gx-1);
    convy = skirtlayers[l].extent.k / (float) (gy-1);

    // if(hght > 0.0f)
    //    cerr << "at " << x << " , " << y << " h = " << hght << endl;
    pos = vpPoint(skirtlayers[l].start.x + (float) x * convx, skirtlayers[l].start.y + hght, skirtlayers[l].start.z + (float) y * convy);*/
}

void Skirts::loadSkirts(string &basefile, int numskirts, int interval)
{
    ifstream skirtfile;
    string skirtfilename;
    float vx, vy, vz, thickness;
    int n_zero = 3;
  
    clear();
    numlayers = numskirts;
    
    for(int s = 0; s < numlayers; s++)
    {
        Skirt skt;
        
        for(int f = 0; f < numframes; f++)
        {
            std::vector<Parcel> parcels;
            int nums;
          
            // string filenum = to_string(f*interval);
            // string filenum = to_string(f);
            string filenum = std::string(n_zero - std::min(n_zero, (int) to_string(f).length()), '0') + to_string(f);
            skirtfilename = basefile + "skirts_" + filenum + "_" + to_string(s) + ".txt";
        
            skirtfile.open((char *) skirtfilename.c_str(), std::ifstream::in);
            if(skirtfile.is_open())
            {
                skirtfile >> nums;
                skirtfile >> thickness;
                thickness *= 1000.0f;
                
                for(int s = 0; s < nums; s++)
                {
                    Parcel p;
                    float c;
                    
                    skirtfile >> vx >> vy >> vz >> c;
                    p.position = vpPoint(vz*1000.0f+30000.0f, vy*1000.0f, vx*1000.0f+30000.0f);
                    p.condensation = c;
                    p.vapour = 0.0f;
                    parcels.push_back(p);
                }
               
                skirtfile.close();
                skt.frames.push_back(parcels);
            }
            else
            {
                cerr << "Error skirt::loadSkirts - unable to open file " << skirtfilename << endl;
                return;
            }
        }
        skirtlayers.push_back(skt);
    }
}

void Plume::loadPlume(string &basefile, int start)
{
    ifstream spheres, densities, charges, sphtypes;
    string sphfilename, dnsfilename, chgfilename, typefilename;
  
    clear();
    
    for(int f = start; f < numframes; f++)
    {
        int sphtcount[4] = {0, 0, 0, 0};
        int n_zero = 3, numspheres;
        std::vector<PlumeSphere> slist;
        string filenum = std::string(n_zero - std::min(n_zero, (int) to_string(f).length()), '0') + to_string(f);
        
        sphfilename = basefile + "spheres_" + filenum + ".dat";
        dnsfilename = basefile + "densities_" + filenum + ".dat";
        chgfilename = basefile + "charge_" + filenum + ".dat";
        typefilename = basefile + "type_" + filenum + ".dat";
        
        // std::ifstream spheres(sphfilename, std::ifstream::binary | std::ifstream::in);
        // std::ifstream densities(dnsfilename, std::ifstream::binary | std::ifstream::in);
        
        // cerr << filename << endl;
        spheres.open((char *) sphfilename.c_str(), std::ifstream::binary | std::ifstream::in);
        if(spheres.is_open())
        {
            // int fsize = (int) spheres.tellg();
            // char * buffer = new char[fsize];
            // spheres.read(buffer, fsize);
            spheres.read((char *) &numspheres, sizeof(int));
            // spheres >> numspheres;
            
            cerr << "Num Spheres = " << numspheres << endl;
            float minr = 100000.0f, maxr = 0.0f;
            float minx = 100000.0f, maxx = -100000.0f, miny = 100000.0f, maxy = -100000.0f, minz = 100000.0f, maxz = -100000.0f;
            
            for(int s = 0; s < numspheres; s++)
            {
                PlumeSphere sph;
                float x, y, z, swp, r;
                
                // spheres >> x >> y >> z >> r >> a;
                spheres.read((char *) &x, sizeof(float));
                spheres.read((char *) &y, sizeof(float));
                spheres.read((char *) &z, sizeof(float));
                spheres.read((char *) &r, sizeof(float));
                //cerr << "[" << f << ", " << s << "] " << x << " " << y << " " << z << " density = " << a << endl;
                x = x * 1000.0f + 30000.0f; // + 15000.0f - 37.0f;
                y = y * 1000.0f; // - 2500.0f + 25.0f;
                z = z * 1000.0f + 30000.0f; // + 15000.0f - 68.0f;
                swp = x; x = z; z = swp;
                sph.pos = vpPoint(x, y, z);
                sph.radius = r * 1000.0f;
                if(x > maxx)
                    maxx = x;
                if(y > maxy)
                    maxy = y;
                if(z > maxz)
                    maxz = z;
                if(x < minx)
                    minx = x;
                if(y < miny)
                    miny = y;
                if(z < minz)
                    minz = z;
                
                if(sph.radius > maxr)
                    maxr = sph.radius;
                if(sph.radius < minr)
                    minr = sph.radius;
                sph.density = 0;
               
                slist.push_back(sph);
            }
            std::cerr << "maxr = " << maxr << " minr = " << minr << std::endl;
            std::cerr << "maxx = " << maxx << " minx = " << minx << std::endl;
            std::cerr << "maxy = " << maxy << " miny = " << miny << std::endl;
            std::cerr << "maxz = " << maxz << " minz = " << minz << std::endl;
            spheres.close();
        }
        else
        {
            cerr << "Error Plume::loadPlume - unable to open file " << sphfilename << endl;
            return;
        }
        
        // cerr << filename << endl;
        densities.open((char *) dnsfilename.c_str(), std::ifstream::binary | std::ifstream::in);
        if(densities.is_open())
        {
            for(int s = 0; s < numspheres; s++)
            {
                float a;
                
                // spheres >> x >> y >> z >> r >> a;
                densities.read((char *) &a, sizeof(float));
                //cerr << "[" << f << ", " << s << "] " << x << " " << y << " " << z << " density = " << a << endl;
                slist[s].density = a;
            }
            densities.close();
        }
        else
        {
            cerr << "Error Plume::loadPlume - unable to open file " << dnsfilename << endl;
            return;
        }
        
        // cerr << filename << endl;
        charges.open((char *) chgfilename.c_str(), std::ifstream::binary | std::ifstream::in);
        if(charges.is_open())
        {
            for(int s = 0; s < numspheres; s++)
            {
                float c;
                
                charges.read((char *) &c, sizeof(float));
                //cerr << "[" << f << ", " << s << "] " << x << " " << y << " " << z << " density = " << a << endl;
                slist[s].charge = c;
                
            }
            charges.close();
        }
        else
        {
            cerr << "Error Plume::loadPlume - unable to open file " << chgfilename << endl;
            return;
        }
        
        sphtypes.open((char *) typefilename.c_str(), std::ifstream::binary | std::ifstream::in);
        if(sphtypes.is_open())
        {
            for(int s = 0; s < numspheres; s++)
            {
                int i;
                
                sphtypes.read((char *) &i, sizeof(int));
                //cerr << "[" << f << ", " << s << "] " << x << " " << y << " " << z << " density = " << a << endl;
                slist[s].sphtype = i;
                sphtcount[i]++;
                
            }
            sphtypes.close();
        }
        else
        {
            cerr << "Error Plume::loadPlume - unable to open file " << typefilename << endl;
            return;
        }
        
        plumespheres.push_back(slist);
        cerr << "Sphere type proportions = " << (float) sphtcount[0] / (float) numspheres << " " <<  (float) sphtcount[1] / (float) numspheres << " " <<  (float) sphtcount[2] / (float) numspheres << " " <<  (float) sphtcount[3] / (float) numspheres << endl;
    }
}
