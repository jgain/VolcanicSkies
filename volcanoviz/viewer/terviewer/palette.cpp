// constraint.cpp: various user generated constraints for ultimate terrain synthesis
// author: James Gain
// date: 5 November 2013
//       21 January 2013 - curve constraints

#include <cassert>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "timer.h"
#include "palette.h"

#include <QtWidgets>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "glwidget.h"

/*
GLfloat manipCol[] = {0.325f, 0.235f, 1.0f, 1.0f};
//GLfloat manipCol[] = {0.406f, 0.294f, 1.0f, 1.0f};
GLfloat curveCol[] = {0.243f, 0.176f, 0.75f, 1.0f};
GLfloat blockedCol[] = {0.5f, 0.5f, 0.8f, 1.0f};
//GLfloat blockedCol[] = {0.325f, 0.235f, 1.0f, 1.0f};
*/

using namespace std;

//
// Palette
//

BrushPalette::BrushPalette(TypeMap * typemap, QWidget *parent)
    : QWidget(parent)
{
    glparent = (GLWidget *) parent;
    int r, c, p = 0;
  
    setAttribute(Qt::WA_StaticContents);
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    tmap = typemap;

    cerr << "2" << endl;
    
    if(!unassignedImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/unassignedIcon.png"))
        cerr << "/unassignedIcon.png not found" << endl;

    if(!activeImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/activeIcon.png"))
        cerr << "/activeIcon.png not found" << endl;
    
    if(!forestImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/forestIcon.png"))
        cerr << "/forestIcon.png not found" << endl;
    
    if(!shrubImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/shrubIcon.png"))
        cerr << "/shrubIcon.png not found" << endl;
    
    if(!obstacleImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/obstacleIcon.png"))
        cerr << "/obstacleIcon.png not found" << endl;

    if(!waterImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/waterIcon.png"))
        cerr << "/waterIcon.png not found" << endl;

    if(!grassImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/grassIcon.png"))
        cerr << "/grassIcon.png not found" << endl;
 
    if(!iceImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/iceIcon.png"))
        cerr <<  "/iceIcon.png not found" << endl;

    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->setColumnStretch(0, 0);
    mainLayout->setColumnStretch(1, 0);
    // mainLayout->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // create as many colour buttons as needed up to limit of PALETTE_ENTRIES
    // Create buttons that are specific to the species
    for(int i = 0; i < PALETTE_ENTRIES; i++)
    {
        if(inSpeciesPalette(i))
        {
            selector[i] = new QPushButton(this);
            r = p / 2; c = p % 2;
            p++;
            currSel = i;
            setDrawType((BrushType) i);
               
            selector[i]->setIconSize(activeImg.size()*1.25f);
            selector[i]->setFixedSize(activeImg.size()*1.25f);
            selector[i]->setFocusPolicy(Qt::NoFocus);
            
            connect(selector[i], &QPushButton::clicked, this, &BrushPalette::typeSelect);
            mainLayout->addWidget(selector[i], r, c, Qt::AlignCenter);
           
        }
    }
    currSel = (int) firstSpeciesPalette();
    setActivePalette();
    setLayout(mainLayout);
    setFocusPolicy(Qt::NoFocus);
}

QSize BrushPalette::sizeHint() const
{
    return QSize(80, 400);
}

void BrushPalette::setDrawType(BrushType btype)
{
    GLfloat * col;
    int r, g, b;
    QString qss;

    typeSel[currSel] = btype;

    cerr << "X" << endl;
    // set colour
    col = tmap->getColour((int) btype);
    cerr << "Y" << endl;
    r = (int) (col[0] * 255.0f);
    g = (int) (col[1] * 255.0f);
    b = (int) (col[2] * 255.0f);
    qss = QString("* { background-color: rgb(%1,%2,%3) }").arg(r).arg(g).arg(b);
    selector[currSel]->setStyleSheet(qss);
    selector[currSel]->show();
    cerr << "Z" << endl;
}

void BrushPalette::setDisabled(BrushType btype)
{
    int r, g, b, bt;
    QString qss;

    bt = (int) btype;

    // set colour to grey
    r = 125; g = 125; b = 125;
    qss = QString("* { background-color: rgb(%1,%2,%3) }").arg(r).arg(g).arg(b);
    selector[bt]->setStyleSheet(qss);
    selector[bt]->show();
    selector[bt]->setEnabled(false);
}

void BrushPalette::setEnabled(BrushType btype)
{
    GLfloat * col;
    int r, g, b, bt;
    QString qss;

    bt = (int) btype;

    // set colour
    col = tmap->getColour(bt);
    r = (int) (col[0] * 255.0f);
    g = (int) (col[1] * 255.0f);
    b = (int) (col[2] * 255.0f);
    qss = QString("* { background-color: rgb(%1,%2,%3) }").arg(r).arg(g).arg(b);
    selector[bt]->setStyleSheet(qss);
    selector[bt]->show();
    selector[bt]->setEnabled(true);
}

BrushType BrushPalette::firstSpeciesPalette()
{
    BrushType fsel;
    fsel = BrushType::FREE;
    return fsel;
}


bool BrushPalette::inSpeciesPalette(int entry)
{
    bool inp = true;
    return inp;
}

void BrushPalette::switchMode(ControlMode cmode)
{
    for(int i = 0; i < PALETTE_ENTRIES; i++)
        if(inSpeciesPalette(i))
            setEnabled((BrushType) i);
    
    /*
    switch(cmode)
    {
        case ControlMode::PAINTBOIDS: // palette disabled
            for(int i = 0; i < PALETTE_ENTRIES; i++)
                if(i < (int) BrushType::BOIDS1)
                    setDisabled((BrushType) i);
                else
                    setEnabled((BrushType) i);
            break;
        case ControlMode::PAINTTYPE: // category palette displayed
            for(int i = 0; i < PALETTE_ENTRIES; i++)
                if(i >= (int) BrushType::GENTLE)
                    setDisabled((BrushType) i);
                else
                    setEnabled((BrushType) i);
            break;
        case ControlMode::PAINTSLOPE: // slope palette displayed
            for(int i = 0; i < PALETTE_ENTRIES; i++)
                if(i < (int) BrushType::GENTLE || i >= (int) BrushType::BOIDS1)
                    setDisabled((BrushType) i);
                else
                    setEnabled((BrushType) i);
            break;
    }*/
}

void BrushPalette::setActivePalette()
{
      
    for(int i = 0; i < PALETTE_ENTRIES; i++)
    {
        if(inSpeciesPalette(i))
        {
            if(i == currSel)
            {
                selector[currSel]->setIcon(QPixmap::fromImage(activeImg));
            }
            else
            {
                selector[i]->setIcon(QIcon());
                if(i == 3)
                    selector[i]->setIcon(QPixmap::fromImage(grassImg));
                if(i == 4)
                    selector[i]->setIcon(QPixmap::fromImage(shrubImg));
                if(i == 5)
                     selector[i]->setIcon(QPixmap::fromImage(forestImg));
                if(i == 6)
                     selector[i]->setIcon(QPixmap::fromImage(waterImg));
                if(i == 8)
                     selector[i]->setIcon(QPixmap::fromImage(iceImg));
                if(i == 9)
                     selector[i]->setIcon(QPixmap::fromImage(obstacleImg));
            }
        }
    }
}

void BrushPalette::directSelect(BrushType brushtype)
{
    currSel = (int) brushtype;
    setActivePalette();
}

void BrushPalette::typeSelect()
{
    for(int i = 0; i < PALETTE_ENTRIES; i++)
    {
        if(inSpeciesPalette(i))
            if(sender() == selector[i])
            {
                currSel = i;
            
                // activate the correct painting mode
                if(i < (int) BrushType::GENTLE)
                    glparent->setMode(ControlMode::PAINTTYPE);
                if(i >= (int) BrushType::GENTLE && i <= (int) BrushType::CLIFF)
                    glparent->setMode(ControlMode::PAINTSLOPE);
                if(i >= (int) BrushType::BOIDS1 && i <= (int) BrushType::BOIDS6)
                    glparent->setMode(ControlMode::PAINTBOIDS);
                if(i >= (int) BrushType::HERDADD && i <= (int) BrushType::HERDDEL)
                    glparent->setMode(ControlMode::PLACEHERD);
            }
    }
    setActivePalette();
}

//
// BrushCursor
//

void BrushCursor::genBrushRing(View * view, Terrain * terrain, float brushradius, bool dashed)
{
    uts::vector<vpPoint> ring;
    int steps, j;
    float a, stepa, tol, tx, ty;
    vpPoint pnt;

   //  shape.clear();
    terrain->getTerrainDim(tx, ty);
    tol = 0.001f * std::max(tx, ty);

    // draw ring to indicate extent of brush stroke
    // generate vertices for ring and drop onto terrain
    if(active)
    {
        a = 0.0f;
        steps = 1000;
        stepa = PI2 / (float) steps;

        for(j = 0; j < steps+1; j++)
        {
            pnt.x = pos.x + cosf(a) * brushradius;
            if(pnt.x >= tx-tolzero) pnt.x = tx-tolzero;
            if(pnt.x <= tolzero) pnt.x = tolzero;
            pnt.y = 0.0f;
            pnt.z = pos.z + sinf(a) * brushradius;
            if(pnt.z >= ty-tolzero) pnt.z = ty-tolzero;
            if(pnt.z <= tolzero) pnt.z = tolzero;
            ring.push_back(pnt);
            a += stepa;
        }
        drapeProject(&ring, &ring, terrain);

        // add height offset to all ring positions
        for(j = 0; j < (int) ring.size(); j++)
            ring[j].y += hghtoffset;

        if(dashed)
            shape.genDashedCylinderCurve(ring, manipradius * 0.5f * view->getScaleFactor(), tol, manipradius * view->getScaleFactor(), 10);
        else
            shape.genCylinderCurve(ring, manipradius * 0.5f * view->getScaleFactor(), tol, 10);

        /*
        // unit sphere at trackpnt
        GLfloat basecol[4];
        basecol[0] = 0.0f;  basecol[1] = 1.0f;  basecol[2] = 0.0f;  basecol[3] = 1.0f;

        debugshape.setColour(basecol);
        glm::mat4 idt = glm::mat4(1.0f), tfm;
        glm::vec3 trs;

        idt = glm::mat4(1.0f);
        trs = glm::vec3(trackpnt.x, trackpnt.y, trackpnt.z);
        tfm = glm::translate(idt, trs);

        debugshape.genSphere(10.0f, 20, 20, tfm);
        */
    }
}

bool BrushCursor::cursorUpdate(View * view, Terrain * terrain, int x, int y, int port[4])
{
    vpPoint frompnt, topnt;

    view->projectingPointPort(x, y, frompnt, port);
    if(terrainProject(frompnt, topnt, view, terrain))
    {
        pos = topnt;
        // trackpnt = pos;
        active = true;
    }
    else
    {
        active = false;
    }
    return active;
}

/// getters and setters for brush radii
void BrushCursor::setRadius(float rad)
{
    radius = rad;
}

//
// BrushPaint
//

BrushPaint::BrushPaint(Terrain * ter, BrushType btype)
{
    terrain = ter;
    brushtype = btype;
    drawing = false;
}

void BrushPaint::paintMap(TypeMap * pmap, float radius)
{
    int dx, dy, si, sj, ei, ej;
    float inr, h, ox, oy, rad;
    vpPoint p;

    terrain->getGridDim(dx, dy);
    pmap->matchDim(dx, dy);

    // apply stroke to type map by setting index values out to a certain radius around the stroke
    rad = terrain->toGrid(radius);
    inr = rad * rad;

    // convert to grid coordinates
    terrain->toGrid(currpnt, ox, oy, h);

    // bound by edge of map
    si = (int) (ox - rad); if(si < 0) si = 0;
    ei = (int) (ox + rad + 0.5f); if(ei >= dx) ei = dx-1;
    sj = (int) (oy - rad); if(sj < 0) sj = 0;
    ej = (int) (oy + rad + 0.5f); if(ej >= dy) ej = dy-1;

    for(int j = sj; j <= ej; j++)
        for(int i = si; i <= ei; i++)
        {
            float cx, cy;

            cx = ox - (float) i; cx *= cx;
            cy = oy - (float) j; cy *= cy;

            if(cx + cy <= inr) // inner region
                pmap->getMap()->set(i, j, (int) brushtype);
        }
}

void BrushPaint::addMousePnt(View * view, TypeMap * pmap, int x, int y, float radius, int port[4])
{
    Region reg;
    bool valid;
    vpPoint prjpnt, terpnt;
    int dx, dy;

    // Timer t;
    // t.start();

    // capture mouse point projected onto terrain
    // must capture current and previous point for cases where updates are not immediate and the mouse has travelled some distance

    view->projectingPointPort(x, y, prjpnt, port);
    // view->projectingPoint(x, y, prjpnt);
    valid = terrainProject(prjpnt, terpnt, view, terrain);

    if(valid)
    {
        if(!drawing) // first point in the stroke
        {
            prevpnt = terpnt;
            currpnt = terpnt;
            drawing = true;
        }
        else
        {
            prevpnt = currpnt;
            currpnt = terpnt;
        }

        // set bounding region to surround the points and their offset radius
        // bnd.reset();
        bnd.includePnt(currpnt);
        bnd.includePnt(prevpnt);

        BoundRect locbnd;
        locbnd = bnd;
        locbnd.expand(radius);
        terrain->getGridDim(dx, dy);

        // convert to terrain coordinates
        reg.x0 = (int) terrain->toGrid(locbnd.min.x);
        if(reg.x0 < 0) reg.x0 = 0;
        reg.y0 = (int) terrain->toGrid(locbnd.min.z);
        if(reg.y0 < 0) reg.y0 = 0;
        reg.x1 = (int) terrain->toGrid(locbnd.max.x);
        if(reg.x1 > dx) reg.x1 = dx;
        reg.y1 = (int) terrain->toGrid(locbnd.max.z);
        if(reg.y1 > dy) reg.y1 = dy;

        paintMap(pmap, radius); // render to paint map
        pmap->setRegion(reg);

        // t.stop();
        // cerr << "Brush time = " << t.peek() << endl;
    }
}

void BrushPaint::startStroke()
{
    bnd.reset();
}

void BrushPaint::finStroke()
{
    drawing = false;
}
