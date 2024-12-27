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


//
// TypeMap
//

#include "typemap.h"
#include "vecpnt.h"
#include "data_importer/data_importer.h"
#include "common/basic_types.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <memory>
#include <functional>

#include <QFileInfo>
#include <QLabel>

/*
Perceptually uniform colourmaps from:
http://peterkovesi.com/projects/colourmaps/
*/

using namespace std;

// palette colours
/*
float freecol[] = {0.755f, 0.645f, 0.538f, 1.0f};
float rockcol[] = {0.2f, 0.2f, 0.2f, 1.0f};
float trailcol[] = {0.882f, 0.843f, 0.713f, 1.0f};
float grasscol[] = {0.749f, 0.815f, 0.611f, 1.0f};
float shrubscol[] = {0.509f, 0.67f, 0.584f, 1.0f};
float forestcol[] = {0.35f, 0.515f, 0.25f, 1.0f};
float watercol[] = {0.537f, 0.623f, 0.752f, 1.0f};
float snowcol[] = {1.0f, 1.0f, 1.0f, 1.0f};
float icecol[] = {0.9f, 0.9f, 1.0f, 1.0f};
float obstaclecol[] = {0.8f, 0.3f, 0.3f, 1.0f};
*/

TypeMap::TypeMap(TypeMapType purpose, int fillval)
{
    tmap = new MapInt;
    setPurpose(purpose);
    if(fillval != -1)
        fill(fillval);
}

TypeMap::TypeMap(int w, int h, TypeMapType purpose, int fillval)
{
    tmap = new MapInt;
    matchDim(w, h);
    setPurpose(purpose);
    if(fillval != -1)
        fill(fillval);
}

TypeMap::~TypeMap()
{
    delete tmap;
    for(int i = 0; i < (int) colmap.size(); i++)
        delete [] colmap[i];
    colmap.clear();
}

void TypeMap::clear()
{
    tmap->fill(0);
}

void TypeMap::initPaletteColTable()
{
    GLfloat *col;

    colmap.clear();
    for(int i = 0; i < 32; i++) // set all colours in table to black initially
    {
        col = new GLfloat[4];
        col[0] = col[1] = col[2] = 0.0f; col[3] = 1.0f;
        colmap.push_back(col);
    }

    numSamples = 21;

    for(int i = 0; i < 4; i++)
    {
        colmap[0][i] = freecol[i];
        colmap[1][i] = icecol[i];
        colmap[2][i] = watercol[i];
        colmap[3][i] = grasscol[i];
        colmap[4][i] = shrubscol[i];
        colmap[5][i] = forestcol[i];
        colmap[6][i] = rockcol[i];
        colmap[7][i] = snowcol[i];
        colmap[8][i] = trailcol[i];
        colmap[9][i] = obstaclecol[i];
        colmap[10][i] =  gentlecol[i];
        colmap[11][i] = steepcol[i];
        colmap[12][i] = cliffcol[i];
      
        colmap[13][i] = boids1col[i];
        colmap[14][i] = boids2col[i];
        colmap[15][i] = boids3col[i];
        colmap[16][i] = boids4col[i];
        colmap[17][i] = boids5col[i];
        colmap[18][i] = boids6col[i];
      
        colmap[19][i] = addcol[i];
        colmap[20][i] = delcol[i];
    }
}

int TypeMap::getNumSamples()
{
    return numSamples;
}

void TypeMap::initPerceptualColTable(std::string colmapfile, int samples, float truncend)
{
    GLfloat *col;
    float r[256], g[256], b[256];
    ifstream infile;
    string valstr, line;
    int i, pos, step;

    if(samples < 3 || samples > 32)
        cerr << "Error: sampling of colour map must be in the range [3,32]" << endl;

    for(i = 0; i < 32; i++) // set all colours in table to black initially
    {
        col = new GLfloat[4];
        col[0] = col[1] = col[2] = 0.0f; col[3] = 1.0f;
        colmap.push_back(col);
    }

    // input is a csv file, with 256 RGB entries, one on each line
    // note that this is not robust to format errors in the input file
    infile.open((char *) colmapfile.c_str(), ios_base::in);

    if(infile.is_open())
    {
        i = 0;
        while(std::getline(infile, line))
        {
            std::size_t prev = 0, pos;

            // red component
            pos = line.find_first_of(",", prev);
            valstr = line.substr(prev, pos-prev);
            istringstream isr(valstr);
            isr >> r[i];
            prev = pos+1;

            // green component
            pos = line.find_first_of(",", prev);
            valstr = line.substr(prev, pos-prev);
            istringstream isg(valstr);
            isg >> g[i];
            prev = pos+1;

            // blue component
            valstr = line.substr(prev, std::string::npos);
            istringstream isb(valstr);
            isb >> b[i];

            i++;
        }
        infile.close();
        numSamples = samples+1;
    }
    else
    {
        std::cerr << "Warning: could not read colourmap at " << colmapfile << std::endl;
        numSamples = -1;
    }

    // now sample the colour map at even intervals according to the number of samples
    // first and last samples map to the beginning and end of the scale
    step = (int) ((256.0f * truncend) / (float) (samples-1));
    pos = 0;
    for(i = 1; i <= samples; i++)
    {
        colmap[i][0] = (GLfloat) r[pos]; colmap[i][1] = (GLfloat) g[pos]; colmap[i][2] = (GLfloat) b[pos];
        pos += step;
    }
}

void TypeMap::clipRegion(Region &reg)
{
    if(reg.x0 < 0) reg.x0 = 0;
    if(reg.y0 < 0) reg.y0 = 0;
    if(reg.x1 > width()) reg.x1 = width();
    if(reg.y1 > height()) reg.y1 = height();
}

void TypeMap::matchDim(int w, int h)
{
    int mx, my;

    mx = tmap->width();
    my = tmap->height();

    // if dimensions don't match then reallocate
    if(w != mx || h != my)
    {
        dirtyreg = Region(0, 0, w, h);
        tmap->setDim(w, h);
        tmap->fill(0); // set to empty type
    }
}

void TypeMap::replaceMap(MapInt * newmap)
{
    assert(tmap->width() == newmap->width());
    assert(tmap->height() == newmap->height());
    for (int y = 0; y < tmap->height(); y++)
        for (int x = 0; x < tmap->width(); x++)
            tmap->set(y,x, newmap->get(y,x));
}

int TypeMap::load(const std::string &filename, TypeMapType purpose)
{
    int tp, maxtp = 0; // mintp = 100;
    int width, height;
    ifstream infile;

    infile.open((char *) filename.c_str(), ios_base::in);
    if(infile.is_open())
    {
        infile >> width >> height;
        // cerr << "width = " << width << " height = " << height << endl;
        matchDim(width, height);
        // convert to internal type map format

        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                switch(purpose)
                {
                    case TypeMapType::BOIDS: // do nothing
                        break;
                    case TypeMapType::CATEGORY:
                    case TypeMapType::SLOPE:
                        infile >> tp;
                        tp++;
                        break;
                   /*
                    case TypeMapType::SUNLIGHT:
                        infile >> val;
                        if(val > maxval)
                            maxval = val;

                        // discretise into ranges of illumination values
                        range = 12.0f; // hours of sunlight
                        // clamp values to range
                        if(val < 0.0f) val = 0.0f;
                        if(val > range) val = range;
                        tp = (int) (val / (range+pluszero) * (numSamples-1))+1;
                        break;*/
                    default:
                        break;
                }
                tmap->set(y,x,tp);

                if(tp > maxtp)
                    maxtp = tp;
            }
        }
        infile.close();
        // cerr << "maxtp = " << maxtp << endl;
        // cerr << "mintp = " << mintp << endl;
    }
    else
    {
        cerr << "Error TypeMap::loadTxt: unable to open file" << filename << endl;
    }
    return maxtp;
}

bool TypeMap::loadCategoryImage(const std::string &filename)
{
    int width, height;
    QImage img(QString::fromStdString(filename)); // load image from file

    QFileInfo check_file(QString::fromStdString(filename));

    if(!(check_file.exists() && check_file.isFile()))
        return false;

    // set internal storage dimensions
    width = img.width();
    height = img.height();
    matchDim(width, height);

    // convert to internal type map format
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++)
        {
            QColor col = img.pixelColor(x, y);
            int r, g, b;
            col.getRgb(&r, &g, &b); // all channels store the same info so just use red
            tmap->set(y,x, r - 100); // convert greyscale colour to category index
        }
    return true;
}

bool TypeMap::readCategoryFile(const std::string &filename)
{
    int val;
    ifstream infile;
    int gx, gy;
    int cnt1=0, cnt2=0, cnt3=0, cnt4=0, cnt5=0, cnt6=0, cnt7=0;

    infile.open((char *) filename.c_str(), ios_base::in);
    if(infile.is_open())
    {
        infile >> gx >> gy;

        if((gx != width()) || (gy != height()))
        {
            cerr << "Error TypeMap::readCategoryFile: map dimensions do not match terrain" << endl;
            cerr << "Terrain = " << width() << " X " << height() << endl;
            cerr << "File = " << gx << " X " << gy << endl;
        }

        for (int y = 0; y < gy; y++)
            for (int x = 0; x < gx; x++)
             {
                infile >> val;
                if(val == 1)
                    cnt1++;
                if(val == 2)
                    cnt2++;
                if(val == 3)
                    cnt3++;
                if(val == 4)
                    cnt4++;
                if(val == 5)
                    cnt5++;
                if(val == 6)
                    cnt6++;
                if(val == 7)
                    cnt7++;
                set(x, y, val);
             }
        infile.close();
        // cerr << "CATEGORY COUNTS: type1 = " << cnt1 << " type2 = " << cnt2 << " type3 = " << cnt3 << endl;
        // cerr << "type4 = " << cnt4 << " type5 = " << cnt5 << " type6 = " << cnt6 << " type7 = " << cnt7 << endl;
        return true;
    }
    else
    {
        cerr << "Error TypeMap::readCategoryFile: unable to open file" << filename << endl;
        return false;
    }
}

int TypeMap::convert(MapFloat * map, TypeMapType purpose, float range)
{
    int tp, maxtp = 0;
    int width, height;

    map->getDim(width, height);
    matchDim(width, height);
    // convert to internal type map format

    for(int x = 0; x < width; x++)
        for(int y = 0; y < height; y++)
        {
            tp = 0;
            
            if(purpose == TypeMapType::SLOPE)
            {
               
                float val = map->get(x, y);
                // tp = (int) (val / (range+pluszero) * (numSamples-1))+1;
                
                
                // discretise into 3 classes according to slope thresholds
                if(val < def_gentleslope)
                {
                    tp = (int) BrushType::GENTLE;
                }
                else
                {
                    if(val < def_steepslope)
                        tp = (int) BrushType::STEEP;
                    else
                        tp = (int) BrushType::CLIFF;
                }
            }
            tmap->set(x,y,tp);

            if(tp > maxtp)
                maxtp = tp;
        }
    /*
    if(purpose == TypeMapType::CDM)
    {
        cerr << "Minimum colour value = " << mincm << endl;
        cerr << "Maxiumum colour value = " << maxcm << endl;
    }*/
    return maxtp;
}

void TypeMap::save(const std::string &filename)
{
    ofstream outfile;

    outfile.open((char *) filename.c_str(), ios_base::out);
    if(outfile.is_open())
    {
        outfile << width() << " " << height() << endl;

        // dimensions
        for (int x = 0; x < width(); x++)
            for (int y = 0; y < height(); y++)
            {
                outfile << get(x, y) << " ";
            }
        outfile.close();
    }
    else
    {
        cerr << "Error TypeMap::save: unable to write to file" << endl;
    }
}

void TypeMap::saveToPaintImage(const std::string &filename)
{
    unsigned char * mask = new unsigned char[tmap->width()*tmap->height()];
    int i = 0;

    cerr << "paint file: " << filename << endl;

    //mask.resize(tmap->width()*tmap->height(), 0.0f);
    for (int x = 0; x < tmap->width(); x++)
        for (int y = 0; y < tmap->height(); y++)
        {
            switch(tmap->get(x,y)) // check order
            {
            case 0:
                mask[i] = 0;
                break;
            case 1: // sparse low
                mask[i] = 38;
                break;
            case 2: // sparse med
                mask[i] = 76;
                break;
            case 3: // sparse tall
                mask[i] = 115;
                break;
            case 4: // dense low
                mask[i] = 153;
                break;
            case 5: // dense med
                mask[i] = 191;
                break;
            case 6: // dense tall
                mask[i] = 230;
                break;
            default:
                mask[i] = 0;
            }
            i++;
        }

    // use QT image save functions
    QImage img;
    img = QImage(mask, tmap->width(), tmap->height(), QImage::Format_Grayscale8);
    img.save(QString::fromStdString(filename), "PNG", 100);
    delete [] mask;
}

void TypeMap::setPurpose(TypeMapType purpose)
{
    usage = purpose;
    switch(usage)
    {
        case TypeMapType::DISPLAY:
            initPaletteColTable();
            break;
        case TypeMapType::BOIDS:
        case TypeMapType::CATEGORY:
        case TypeMapType::SLOPE:
            /*
            initPerceptualColTable("../../common/colourmaps/isoluminant_cm_70_c39_n256.csv", 20);
            // replace 0 with natural terrain colour
            colmap[1][0] = 0.7f; colmap[1][1] = 0.6f; colmap[1][2] = 0.5f; // transparent
            colmap[2][0] = 0.0f; colmap[2][1] = 0.0f; colmap[2][2] = 1.0f; // black
            colmap[numSamples+2][0] = 1.0f; colmap[numSamples+2][1] = 0.0f; colmap[numSamples+2][2] = 0.0f; // red
             */
            break;
        case TypeMapType::TRAIL:
        case TypeMapType::HOME:
            initPaletteColTable();
            break;
        default:
            cerr << "Purpose not correctly specified" << endl;
            break;
    }
}

void TypeMap::resetType(int ind)
{
    // wipe all previous occurrences of ind
    #pragma omp parallel for
    for(int j = 0; j < tmap->height(); j++)
        for(int i = 0; i < tmap->width(); i++)
            if(tmap->get(j,i) == ind)
                tmap->set(j,i,0);
    dirtyreg.x0 = 0; dirtyreg.y0 = 0;
    dirtyreg.x1 = tmap->width(); dirtyreg.y1 = tmap->height();
}

void TypeMap::setColour(int ind, GLfloat * col)
{
    for(int i = 0; i < 4; i++)
        colmap[ind][i] = col[i];
}
