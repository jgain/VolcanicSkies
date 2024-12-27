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


// terrain.h: model for terrain. Responsible for storage and display of heightfield terrain data
// author: James Gain
// date: 17 December 2012

#ifdef _WIN32
#include <glew.h>
#else
#include <GL/glew.h>
#endif

#include "terrain.h"
#include <sstream>
#include <streambuf>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <sys/stat.h>


using namespace std;

void Terrain::toGrid(vpPoint p, float & x, float & y, float & h) const
{
    int gx, gy;
    float tx, ty, convx, convy;

    getGridDim(gx, gy);
    getTerrainDim(tx, ty);
    convx = (float) (gx-1) / tx;
    convy = (float) (gy-1) / ty;
    x = p.x * convx;
    y = p.z * convy;
    h = p.y;
    if(scaleOn)
        h *= scfac;
}


void Terrain::toGrid(vpPoint p, int &x, int &y) const
{
    int gx, gy;
    float tx, ty, convx, convy;

    getGridDim(gx, gy);
    getTerrainDim(tx, ty);
    convx = (float) (gx-1) / tx;
    convy = (float) (gy-1) / ty;
    x = (int) (p.x * convx);
    y = (int) (p.z * convy);
}


float Terrain::toGrid(float wdist) const
{
    int gx, gy;
    float tx, ty, conv;

    getGridDim(gx, gy);
    getTerrainDim(tx, ty);
    conv = (float) (gx-1) / tx;
    return wdist * conv;
}


vpPoint Terrain::toWorld(float x, float y, float h) const
{
    int gx, gy;
    float tx, ty, convx, convy;

    getGridDim(gx, gy);
    getTerrainDim(tx, ty);
    convx = tx / (float) (gx-1);
    convy = ty / (float) (gy-1);

    return vpPoint(x * convx, h, y * convy);
}

vpPoint Terrain::toWorld(int x, int y, float h) const
{
    int gx, gy;
    float tx, ty, convx, convy;

    getGridDim(gx, gy);
    getTerrainDim(tx, ty);
    convx = tx / (float) (gx-1);
    convy = ty / (float) (gy-1);

    return vpPoint((float) x * convx, h, (float) y * convy);
}


float Terrain::toWorld(float gdist) const
{
    int gx, gy;
    float tx, ty, conv;

    getGridDim(gx, gy);
    getTerrainDim(tx, ty);
    conv = tx / (float) (gx-1);

    return gdist * conv;
}


bool Terrain::inWorldBounds(vpPoint p) const
{
    return (p.x >= 0.0f && p.x <= dimx && p.z >= 0.0f && p.z <= dimy);
}

bool Terrain::inSynthBounds(vpPoint p) const
{
    return (p.x >= 0.0f-synthx && p.x <= dimx+synthx && p.z >= 0.0f-synthy && p.z <= dimy+synthy);
}

void Terrain::init(int dx, int dy, float sx, float sy)
{
    grid->setDim(dx, dy);
    setTerrainDim(sx, sy);
    setFocus(vpPoint(sx/2.0f, grid->get(dy/2-1,dx/2-1), sy/2.0f));
    scfac = 1.0f;

    // init accel structure
    spherestep = 8;
    numspx = (grid->width()-1) / spherestep + 1; numspy = (grid->height()-1) / spherestep + 1;
    for(int i = 0; i < numspx; i++)
    {
        std::vector<AccelSphere> sphrow;
        for(int j = 0; j < numspy; j++)
        {
            AccelSphere sph;
            sphrow.push_back(sph);
        }
        boundspheres.push_back(sphrow);
    }

    bufferState = BufferState::REALLOCATE;
    accelValid = false;
    scaleOn = false;
}

void Terrain::initGrid(int dx, int dy, float sx, float sy)
{
    init(dx, dy, sx, sy);
    grid->fill(0.0f);
}

void Terrain::delGrid()
{
    if(boundspheres.size() > 0)
    {
        for(int i = 0; i < (int) boundspheres.size(); i++)
            boundspheres[i].clear();
        boundspheres.clear();
    }

    bufferState = BufferState::REALLOCATE;
    accelValid = false;
}

void Terrain::clipRegion(Region &reg)
{
    if(reg.x0 < 0) reg.x0 = 0;
    if(reg.y0 < 0) reg.y0 = 0;
    if(reg.x1 > grid->width()) reg.x1 = grid->width();
    if(reg.y1 > grid->height()) reg.y1 = grid->height();
}

void Terrain::setMidFocus()
{
    int dx, dy;
    float sx, sy;
    
    getGridDim(dx, dy);
    getTerrainDim(sx, sy);
    if(dx > 0 && dy > 0)
        setFocus(vpPoint(sx/2.0f, grid->get(dy/2-1,dx/2-1), sy/2.0f));
    else
        setFocus(vpPoint(0.0f, 0.0f, 0.0f));
}

void Terrain::getMidPoint(vpPoint & mid)
{
    int dx, dy;
    float sx, sy;

    getGridDim(dx, dy);
    getTerrainDim(sx, sy);
    if(dx > 0 && dy > 0)
        mid = vpPoint(sx/2.0f, grid->get(dy/2-1,dx/2-1), sy/2.0f);
    else
        mid = vpPoint(0.0f, 0.0f, 0.0f);

}

void Terrain::getGridDim(int & dx, int & dy) const
{
    dx = grid->width();
    dy = grid->height();
}

void Terrain::getGridDim(uint & dx, uint & dy)
{
    dx = (uint) grid->width();
    dy = (uint) grid->height();
}

void Terrain::getTerrainDim(float &tx, float &ty) const
{
    tx = dimx; ty = dimy;
}

void Terrain::setTerrainDim(float tx, float ty)
{
    int gx, gy;

    getGridDim(gx, gy);
    dimx = tx; dimy = ty;

    // calc allowable synth border
    synthx = (0.5f / (float) (gx-1) + pluszero) * dimx;
    synthy = (0.5f / (float) (gy-1) + pluszero) * dimy;
}

float Terrain::getTerrainHectArea()
{
    float tx, ty, area;

    getTerrainDim(tx, ty);
    area = tx * ty; // sq metres
    area /= 10000.0f; // 10,000 sq metres in a hectare
    return area;
}

float Terrain::samplingDist()
{
    int dx, dy;
    float tx, ty;
    getGridDim(dx, dy);
    getTerrainDim(tx, ty);
    return (0.5f * std::min(tx, ty)) / (float) (std::max(dx,dy)-1); // separation between vertices, about 2-3 vertices per grid cell
}

float Terrain::smoothingDist()
{
    return 30.0f * samplingDist(); // about 10 grid points per curve segment
}

float Terrain::longEdgeDist()
{
    float tx, ty;
    getTerrainDim(tx, ty);
    return std::max(tx, ty);
}

float Terrain::getHeight(int x, int y)
{
    return grid->get(x,y);
}

float Terrain::getFlatHeight(int idx)
{
    int x, y, dx, dy;
    getGridDim(dx, dy);

    x = idx % dx;
    y = idx / dx;
    return grid->get(x,y);
}

void Terrain::getNormal(int x, int y, Vector & norm, int reach)
{
    vpPoint x1, x2, y1, y2;
    int dx, dy;
    Vector dfdx, dfdy;

    getGridDim(dx, dy);

    // x-positions
    if(x > reach-1)
        x1 = toWorld(x-reach, y, getHeight(x-reach, y));
    else
        x1 = toWorld(0, y, getHeight(0, y));

    if(x < dx-reach)
        x2 = toWorld(x+reach, y, getHeight(x+reach, y));
    else
        x2 = toWorld(dx-1, y, getHeight(dx-1, y));

    // y-positions
    if(y > reach-1)
        y1 = toWorld(x, y-reach, getHeight(x, y-reach));
    else
        y1 = toWorld(x, 0, getHeight(x, 0));

    if(y < dy-reach)
        y2 = toWorld(x, y+reach, getHeight(x, y+reach));
    else
        y2 = toWorld(x, dy-1, getHeight(x, dy-1));

    // cross pattern
    dfdx.diff(x1, x2);
    dfdy.diff(y1, y2);
    dfdx.normalize();
    dfdy.normalize();

    norm.cross(dfdx, dfdy);
    norm.mult(-1.0f);
    //norm.normalize();
}

float Terrain::getSlope(int x, int y)
{
    Vector norm, up;

    getNormal(x, y, norm, 1);
    up = Vector(0.0f, 1.0f, 0.0f);
    return norm.angle(up);
}

float Terrain::getCellExtent()
{
    return dimx / (float) grid->width();
}

void Terrain::updateBuffers(PMrender::TRenderer * renderer) const
{
    const int width = grid->width();
    const int height = grid->height();
    float scx, scy;

    getTerrainDim(scx, scy);

    glewExperimental = GL_TRUE;
    if(!glewSetupDone)
      {
         GLenum err = glewInit();
         if (GLEW_OK != err)
         {
            std::cerr<< "GLEW: initialization failed\n\n";
         }
         glewSetupDone = true;
      }

    if (bufferState == BufferState::REALLOCATE || bufferState == BufferState::DIRTY )
        renderer->updateHeightMap(width, height, scx, scy, (GLfloat*)grid->getPtr(), true);
    else
        renderer->updateHeightMap(width, height, scx, scy, (GLfloat*)grid->getPtr());

    bufferState = BufferState::CLEAN;
}

void Terrain::draw(View * view, PMrender::TRenderer *renderer) const
{
    updateBuffers(renderer);

    // call draw function
    renderer->draw(view);
}

void Terrain::buildSphereAccel()
{
    int si, sj, i, j, imin, imax, jmin, jmax;
    float rad, sqlen;
    vpPoint p, c, b1, b2, d;
    Vector del;

    // cerr << "numspx = " << numspx << ", numspy = " << numspy << endl;
    for(si = 0; si < numspx; si++)
        for(sj = 0; sj < numspy; sj++)
        {
            imin = si*spherestep; imax = std::min(imin+spherestep, grid->width());
            jmin = sj*spherestep; jmax = std::min(jmin+spherestep, grid->height());
            // cerr << "(" << si << ", " << sj << ") = " << "i: " << imin << " - " << imax << " j: " << jmin << " - " << jmax << endl;

            // center point 
            b1 = toWorld(imin, jmin, 0.0f); // grid->get(jmin,imin));
            d = vpPoint(b1.z, b1.y, b1.x);
            drapePnt(d, d);
            b1.y = d.y;
            b2 = toWorld(imax, jmax, 0.0f); // , grid->get(jmax-1,imax-1));
            d = vpPoint(b2.z, b2.y, b2.x);
            drapePnt(d, d);
            b2.y = d.y;
            c.affinecombine(0.5f, b1, 0.5f, b2);

            // update radius
            rad = 0.0f;
            for(j = jmin; j < jmax; j++)
                for(i = imin; i < imax; i++)
                {
                    p = toWorld(i, j, grid->get(j,i));
                    del.diff(c, p);
                    sqlen = del.sqrdlength();
                    if(sqlen > rad)
                        rad = sqlen;
                }
            boundspheres[si][sj].center = c;
            boundspheres[si][sj].radius = sqrtf(rad);
        }
    accelValid = true;
}

bool Terrain::rayIntersect(vpPoint start, Vector dirn, vpPoint & p)
{
    int i, j, si, sj, imin, imax, jmin, jmax;
    vpPoint currp;
    float besttval, tval, dist;
    bool found = false;
    float tol = dimx / (float) (grid->width()-1); // set world space detection tolerance to approx half gap between grid points
    int dx, dy;

    // tol *= 5.0f;

    getGridDim(dx, dy);
    besttval = 100000000.0f;

    /*
    for(int i = 0; i < dx; i++)
        for(int j = 0; j < dy; j++)
        {
            currp = toWorld(i, j, 0.0f);
            vpPoint dpnt = vpPoint(currp.z, currp.y, currp.x);
            drapePnt(dpnt, dpnt);
            currp.y = dpnt.y;

            rayPointDist(start, dirn, currp, tval, dist);
            if(dist < tol)
            {
                found = true;
                if(tval < besttval)
                {
                    besttval = tval;
                    p = currp;
                }
            }
        }
    */

    if(!accelValid)
      buildSphereAccel();

    // bounding sphere accel structure
    for(si = 0; si < numspx; si++)
        for(sj = 0; sj < numspy; sj++)
        {
            rayPointDist(start, dirn, boundspheres[si][sj].center, tval, dist);
            if(dist <= boundspheres[si][sj].radius) // intersects enclosing sphere so test enclosed points
            {
                imin = si*spherestep; imax = std::min(imin+spherestep, grid->width());
                jmin = sj*spherestep; jmax = std::min(jmin+spherestep, grid->height());
                // check ray against grid points
                for(j = jmin; j < jmax; j++)
                    for(i = imin; i < imax; i++)
                    {
                        currp = toWorld(i, j, 0.0f);
                        vpPoint dpnt = vpPoint(currp.z, currp.y, currp.x);
                        drapePnt(dpnt, dpnt);
                        currp.y = dpnt.y;

                        rayPointDist(start, dirn, currp, tval, dist);
                        if(dist < tol)
                        {
                            found = true;
                            if(tval < besttval)
                            {
                                besttval = tval;
                                p = currp;
                            }
                        }
                    }

            }

        }
    return found;
}

bool Terrain::pick(int sx, int sy, View * view, vpPoint & p)
{
    vpPoint start;
    Vector dirn;

    cerr << "sx = " << sx << ", sy = " << sy << endl;

    // find ray params from viewpoint through screen <sx, sy>
    view->projectingRay(sx, sy, start, dirn);

    return rayIntersect(start, dirn, p);
}

bool Terrain::drapePnt(vpPoint pnt, vpPoint & drape)
{
    float x, y, h, drapeh, u, v, h0, h1, ux, uy;
    int cx, cy, dx, dy;

    getGridDim(dx, dy);
    toGrid(pnt, x, y, h); // locate point on base domain

    // test whether point is in bounds
    ux = (float) (dx-1) - pluszero;
    uy = (float) (dy-1) - pluszero;

    if(x < pluszero || y < pluszero || x > ux || y > uy)
        return false;

    // index of grid cell
    cx = (int) floor(x);
    cy = (int) floor(y);

    // get parametric coordinates within grid cell
    u = (x - (float) cx);
    v = (y - (float) cy);

    // bilinear interpolation
    h0 = (1.0f - u) * grid->get(cy,cx) + u * grid->get(cy,cx+1);
    h1 = (1.0f - u) * grid->get(cy+1,cx) + u * grid->get(cy+1,cx);
    drapeh = (1.0f - v) * h0 + v * h1;
    // this could be implemented using ray-triangle intersection
    // but it would be much less efficient
    drape = toWorld(x, y, drapeh);

    return true;
}

long Terrain::hashVert(vpPoint pnt, BoundBox bbox)
{
    long x, y, z;
    float range = 2500.0f;
    long lrangesq, lrange = 2500;

    lrangesq = lrange * lrange;

    // discretise vertex within bounds of the enclosing bounding box
    x = (long) (((pnt.x - bbox.min.x) * range) / bbox.diagLen()) * lrangesq;
    y = (long) (((pnt.y - bbox.min.y) * range) / bbox.diagLen()) * lrange;
    z = (long) (((pnt.z - bbox.min.z) * range) / bbox.diagLen());
    return x+y+z;
}

bool Terrain::loadSTL(string &filename)
{
    ifstream infile;
    char * inbuffer;
    struct stat results;
    int insize, inpos, numt, t, i;
    vpPoint vpos;
    std::vector<vpPoint> verts;
    // Triangle tri;

    // assumes binary format STL file
    infile.open((char *) filename.c_str(), ios_base::in | ios_base::binary);
    if(infile.is_open())
    {
        // get the size of the file
        stat((char *) filename.c_str(), &results);
        insize = results.st_size;

        // put file contents in buffer
        inbuffer = new char[insize];
        infile.read(inbuffer, insize);
        if(!infile) // failed to read from the file for some reason
        {
            cerr << "Error Mesh::readSTL: unable to populate read buffer" << endl;
            return false;
        }

        // interpret buffer as STL file
        if(insize <= 84)
        {
            cerr << "Error Mesh::readSTL: invalid STL binary file, too small" << endl;
            return false;
        }

        inpos = 80; // skip 80 character header
        if(inpos+4 >= insize){ cerr << "Error Mesh::readSTL: malformed header on stl file" << endl; return false; }
        numt = (int) (* ((long *) &inbuffer[inpos]));
        inpos += 4;

        t = 0;

        // triangle vertices have consistent outward facing clockwise winding (right hand rule)
        while(t < numt) // read in triangle data
        {
            // normal
            if(inpos+12 >= insize){ cerr << "Error Mesh::readSTL: malformed stl file" << endl; return false; }
            // IEEE floating point 4-byte binary numerical representation, IEEE754, little endian
            // tri.n = cgp::Vector((* ((float *) &inbuffer[inpos])), (* ((float *) &inbuffer[inpos+4])), (* ((float *) &inbuffer[inpos+8])));
            inpos += 12;

            // vertices
            for(i = 0; i < 3; i++)
            {
                if(inpos+12 >= insize){ cerr << "Error Mesh::readSTL: malformed stl file" << endl; return false; }
                vpos = vpPoint((* ((float *) &inbuffer[inpos])), (* ((float *) &inbuffer[inpos+4])), (* ((float *) &inbuffer[inpos+8])));
                // tri.v[i] = (int) verts.size();
                verts.push_back(vpos);
                inpos += 12;
            }
            // tris.push_back(tri);
            t++;
            inpos += 2; // handle attribute byte count - which can simply be discarded
        }

        // tidy up
        delete inbuffer;
        infile.close();

        cerr << "num vertices = " << (int) verts.size() << endl;
        // cerr << "num triangles = " << (int) tris.size() << endl;

        // STL provides a triangle soup so merge vertices that are coincident
        vector<vpPoint> cleanverts;
        long key;
        int i, hitcount = 0;
        // use hashmap to quickly look up vertices with the same coordinates
        std::unordered_map<long, int> idxlookup; // key is concatenation of vertex position, value is index into the cleanverts vector
        BoundBox bbox;

        // construct a bounding box enclosing all vertices
        for(i = 0; i < (int) verts.size(); i++)
            bbox.includePnt(verts[i]);

        // remove duplicate vertices
        for(i = 0; i < (int) verts.size(); i++)
        {
            key = hashVert(verts[i], bbox);
            if(idxlookup.find(key) == idxlookup.end()) // key not in map
            {
                idxlookup[key] = (int) cleanverts.size(); // put index in map for quick lookup
                cleanverts.push_back(verts[i]);
            }
            else
            {
                hitcount++;
            }
        }
        cerr << "num duplicate vertices found = " << hitcount << " of " << (int) verts.size() << endl;
        cerr << "clean verts = " << (int) cleanverts.size() << endl;
        cerr << "bbox min = " << bbox.min.x << ", " << bbox.min.y << ", " << bbox.min.z << endl;
        cerr << "bbox max = " << bbox.max.x << ", " << bbox.max.y << ", " << bbox.max.z << endl;
        cerr << "bbox diag = " << bbox.diagLen() << endl;

        verts.clear();
        verts = cleanverts;
        
        for(i = 0; i < 20; i++)
        {
            cerr << verts[i].x << " " << verts[i].y << " " << verts[i].z << endl;
        }
    }
    else
    {
        cerr << "Error Terrain::readSTL: unable to open " << filename << endl;
        return false;
    }
    return true;

}

void Terrain::loadElv(const uts::string &filename, TypeMap * tmap)
{
    float lat;
    int dx, dy;

    float val;
    // int code;
    ifstream infile;

    infile.open((char *) filename.c_str(), ios_base::in);
    if(infile.is_open())
    {
        infile >> dx >> dy;
        infile >> step;
        infile >> lat;
        delGrid();
        init(dx, dy, (float) dx * step, (float) dy * step);
        tmap->matchDim(dx, dy);
        latitude = lat;
        for (int y = 0; y < dy; y++)
        {
            for (int x = 0; x < dx; x++)
            {

                infile >> val;
                grid->set(x,y, val); // (val * 0.3048f)); // convert from feet to metres
            }
        }
        
        /*
        for (int x = 0; x < dx; x++)
        {
            for (int y = 0; y < dy; y++)
            {

                infile >> code;
                tmap->set(x, y, code);
            }
        }*/
        
        /*
        for (int x = 0; x < dx; x++)
        {
            for (int y = 0; y < dy; y++)
            {

                infile >> val;
                
                // vmap->set(x, y, val);
                int numSamples = 20;
                code = (int) (val / (1.0f+pluszero) * (numSamples-1))+1;
                tmap->set(x, y, code);
                // cerr << code << " ";
            }
        }*/
        setMidFocus();
        infile.close();
    }
    else
    {
        cerr << "Error Terrain::loadElv:unable to open file " << filename << endl;
    }
}

// t is a value that goes from 0 to 1 to interpolate in a C1 continuous way across uniformly sampled data points.
// when t is 0, this will return B.  When t is 1, this will return C.  Inbetween values will return an interpolation
// between B and C.  A and B are used to calculate slopes at the edges.
float Terrain::cubicHermite (float A, float B, float C, float D, float t)
{
    float a = -A / 2.0f + (3.0f*B) / 2.0f - (3.0f*C) / 2.0f + D / 2.0f;
    float b = A - (5.0f*B) / 2.0f + 2.0f*C - D / 2.0f;
    float c = -A / 2.0f + C / 2.0f;
    float d = B;

    return a*t*t*t + b*t*t + c*t + d;
}

float Terrain::getClamp(MapFloat * map, int x, int y)
{
    int cx, cy, dx, dy;

    map->getDim(dx, dy); dx -= 1; dy -= 1;

    cx = std::max(0,x);
    cx = std::min(dx, cx);
    cy = std::max(0, y);
    cy = std::min(dy, cy);

    return map->get(cx, cy);
}

float Terrain::sampleBicubic (MapFloat * map, int dx, int dy, float u, float v)
{
    // calculate coordinates -> also need to offset by half a pixel to keep image from shifting down and left half a pixel
    float x = (u * dx) - 0.5;
    int xint = int(x);
    float xfract = x - floor(x);

    float y = (v * dy) - 0.5;
    int yint = int(y);
    float yfract = y - floor(y);

    // 1st row
    float p00 = getClamp(map, xint - 1, yint - 1);
    float p10 = getClamp(map, xint + 0, yint - 1);
    float p20 = getClamp(map, xint + 1, yint - 1);
    float p30 = getClamp(map, xint + 2, yint - 1);

    // 2nd row
    float p01 = getClamp(map, xint - 1, yint + 0);
    float p11 = getClamp(map, xint + 0, yint + 0);
    float p21 = getClamp(map, xint + 1, yint + 0);
    float p31 = getClamp(map, xint + 2, yint + 0);

    // 3rd row
    float p02 = getClamp(map, xint - 1, yint + 1);
    float p12 = getClamp(map, xint + 0, yint + 1);
    float p22 = getClamp(map, xint + 1, yint + 1);
    float p32 = getClamp(map, xint + 2, yint + 1);

    // 4th row
    float p03 = getClamp(map, xint - 1, yint + 2);
    float p13 = getClamp(map, xint + 0, yint + 2);
    float p23 = getClamp(map, xint + 1, yint + 2);
    float p33 = getClamp(map, xint + 2, yint + 2);

    // interpolate bi-cubically!
    float col0 = cubicHermite(p00, p10, p20, p30, xfract);
    float col1 = cubicHermite(p01, p11, p21, p31, xfract);
    float col2 = cubicHermite(p02, p12, p22, p32, xfract);
    float col3 = cubicHermite(p03, p13, p23, p33, xfract);
    float value = cubicHermite(col0, col1, col2, col3, yfract);
    return value;
}

void Terrain::genTestElv(const uts::string &outfilename)
{
    float val;
    ofstream outfile;

    outfile.open((char *) outfilename.c_str(), ios_base::out);

    if(outfile.is_open())
    {
        outfile << "512 512 1.0 45.0" << endl;
        val = 512.0f;
        for (int y = 0; y < 512; y++)
        {
            // val = 512.0f;
            for (int x = 0; x < 512; x++)
            {
                outfile << val << " ";

                if(x < 256)
                    val -= 0.5f / 0.3048f;
                else
                    val += 0.5f / 0.3048f;

            }

            if(y < 256)
                val -= 0.5f / 0.3048f;
            else
                val += 0.5f / 0.3048f;

        }
        outfile.close();
    }
    else
    {
        cerr << "Error Terrain::genTestElv:unable to open file " << outfilename << endl;
    }
}

void Terrain::cropElv(const uts::string &infilename, const uts::string &outfilename, int startx, int starty, int extentx, int extenty, int upsample)
{
    float lat;
    int dx, dy;
    MapFloat * termap = new MapFloat, * cropmap = new MapFloat;
    bool valid;

    float val;
    ifstream infile;
    ofstream outfile;

    infile.open((char *) infilename.c_str(), ios_base::in);
    outfile.open((char *) outfilename.c_str(), ios_base::out);

    if(infile.is_open() && outfile.is_open())
    {
        // read in terrain
        infile >> dx >> dy;
        termap->setDim(dx, dy);
        infile >> step;
        infile >> lat;
        delGrid();
        init(dx, dy, (float) dx * step, (float) dy * step);
        latitude = lat;
        for (int y = 0; y < dy; y++)
        {
            for (int x = 0; x < dx; x++)
            {
                infile >> val;
                termap->set(x,y, val);
            }
        }

        // check crop parameters
        valid = startx >= 0 && startx+extentx < dx && starty >= 0 && starty+extenty < dy;
        if(valid)
        {
            // generate crop map
            cropmap->setDim(extentx, extenty);
            for (int y = 0; y < extenty; y++)
                for (int x = 0; x < extentx; x++)
                    cropmap->set(x, y, termap->get(startx+x, starty+y));

            // upsample and save
            outfile << upsample * extentx << " " << upsample * extenty << " " << step / (float) upsample << " " << latitude << endl;
            for (int y = 0; y < upsample * extenty; y++)
            {
                for (int x = 0; x < upsample * extentx; x++)
                {
                    // determine cell and gather surrounding points
                    float u = (float) x / (float) (upsample*extentx-1);
                    float v = (float) y / (float) (upsample*extenty-1);
                    float val = sampleBicubic(cropmap, extentx, extenty, u, v);
                    outfile << val << " ";
                }
            }
        }
        else
        {
            cerr << "Error Terrain::cropElv: invalid parameters out of bounds of terrain" << endl;
        }

        infile.close(); outfile.close();
    }
    else
    {
        cerr << "Error Terrain::cropElv:unable to open file " << infilename << " or " << outfilename << endl;
    }

}


void Terrain::saveElv(const uts::string &filename)
{
    int gx, gy;
    ofstream outfile;

    outfile.open((char *) filename.c_str(), ios_base::out);
    if(outfile.is_open())
    {
        getGridDim(gx, gy);
        outfile << gx << " " << gy << " " << step << " " << latitude << endl;
        for (int x = 0; x < gx; x++)
        {
            for (int y = 0; y < gy; y++)
            {
                outfile << grid->get(x,y) << " ";
            }
        }
        outfile << endl;
        outfile.close();
    }
    else
    {
        cerr << "Error Terrain::loadElv:unable to open file " << filename << endl;
    }
}

void Terrain::calcMeanHeight()
{
    int i, j, cnt = 0;
    hghtmean = 0.0f;

    for(j = 0; j < grid->height(); j++)
        for(i = 0; i < grid->width(); i++)
        {
            hghtmean += grid->get(j,i);
            cnt++;
        }
    hghtmean /= (float) cnt;
}

void Terrain::getHeightBounds(float &minh, float &maxh)
{
    int i, j;
    float hght;

    maxh = -10000000.0f;
    minh = 100000000.0;

    for(j = 0; j < grid->height(); j++)
        for(i = 0; i < grid->width(); i++)
        {
            hght = grid->get(j,i);
            if(hght < minh)
                minh = hght;
            if(hght > maxh)
                maxh = hght;
        }
}
