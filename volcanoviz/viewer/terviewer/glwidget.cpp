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


/****************************************************************************
**
** Copyright (C) 2012 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "glwidget.h"
#include "window.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <QGridLayout>
#include <QGLFramebufferObject>
#include <QImage>
#include <QCoreApplication>
#include <QMessageBox>
#include <QInputDialog>

#include <fstream>

using namespace std;

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

////
// Scene
////

Scene::Scene()
{
    view = new View();
    terrain = new Terrain();
    airlayers = new AirLayers();
    skirts = new Skirts();
    plume = new Plume();
 
    terrain->initGrid(1024, 1024, 10000.0f, 10000.0f);
    view->setForcedFocus(terrain->getFocus());
    view->setViewScale(terrain->longEdgeDist());

    int dx, dy;
    terrain->getGridDim(dx, dy);

    for(TypeMapType t: all_typemaps)
        maps[(int) t] = new TypeMap(dx, dy, t);
    // maps[2]->setRegion(terrain->coverRegion());
    overlay = TypeMapType::DISPLAY;
}

Scene::~Scene()
{
    delete view;
    delete terrain;
    delete skirts;
    delete plume;
    for(TypeMapType t: all_typemaps)
        if(maps[(int) t] != nullptr)
        {
            delete maps[(int) t];
            maps[(int) t] = nullptr;
        }
}

////
// GLWidget
////

GLWidget::GLWidget(const QGLFormat& format, string datadir, QWidget *parent)
    : QGLWidget(format, parent)
{
    this->datadir = datadir;
 
    qtWhite = QColor::fromCmykF(0.0, 0.0, 0.0, 0.0);
    vizpopup = new QLabel();

    // set up animation timers
    atimer = new QTimer(this);
    connect(atimer, SIGNAL(timeout()), this, SLOT(animUpdate()));
    rtimer = new QTimer(this);
    connect(rtimer, SIGNAL(timeout()), this, SLOT(rotateUpdate()));

    glformat = format;

    // main design scene
    addScene();

    currscene = 0;

    mainparent = (Window *) parent;
    renderer = new PMrender::TRenderer(NULL, "../terviewer/shaders/");
    // palette = new BrushPalette(getTypeMap(TypeMapType::CATEGORY), this);
    // setRadius(250.0f);
    cmode = ControlMode::VIEW;
    
    viewing = false;
    viewlock = false;
    decalsbound = false;
    focuschange = false;
    focusviz = false;
    timeron = false;
    dbloaded = false;
    active = true;
    plumevisibility = true;
    skirtsvisibility = true;
    plumecharge = false;

    scf = 10000.0f;
    decalTexture = 0;

    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    resize(sizeHint());
    setSizePolicy (QSizePolicy::Ignored, QSizePolicy::Ignored);
    
   
    GLfloat tcol[4], scol[4], pcol[4], acol[4];
    
    genAirLayerShape();
    acol[0] = 1.0f; acol[1] = 1.0f; acol[2] = 1.0f; acol[3] = 1.0f;
    airlayershape.setColour(acol);
    
    genSkirtShape();
    scol[0] = 1.0f; scol[1] = 1.0f; scol[2] = 1.0f; scol[3] = 1.0f;
    skirtshape.setColour(scol);
    
    genPlumeShape();
    for(int i = 0; i < 3; i++)
        pcol[i] = 0.2f;
    pcol[3] = 1.0f;
    plumeshape.setColour(pcol);
    
    genTreeShape(0.2, 0.1);
    for(int i = 0; i < 4; i++)
        tcol[i] = forestcol[i];
    treeshape.setColour(tcol);
    
    genShrubShape(0.1, 0.2);
    for(int i = 0; i < 4; i++)
        scol[i] = shrubscol[i];
    shrubshape.setFlatColour(scol);
    
    nfield1 = new NoiseField(getTerrain(), 2, 0);
    nfield2 = new NoiseField(getTerrain(), 3, 1198724);
    
    update();
}

GLWidget::~GLWidget()
{
    delete atimer;
    delete rtimer;
    if(vizpopup) delete vizpopup;
    if (renderer) delete renderer;

    // delete views
    for(int i = 0; i < (int) scenes.size(); i++)
        delete scenes[i];

    if (decalTexture != 0)	glDeleteTextures(1, &decalTexture);
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(100, 80);
}

QSize GLWidget::sizeHint() const
{
    return QSize(1000, 800);
}


void GLWidget::screenCapture(QImage * capImg, QSize capSize)
{
    paintGL();
    glFlush();

    (* capImg) = grabFrameBuffer();
    (* capImg) = capImg->scaled(capSize, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
}

void GLWidget::unscaledScreenCapture(QImage * capImg)
{
    paintGL();
    glFlush();

    (* capImg) = grabFrameBuffer();
}

View * GLWidget::getView()
{
    if((int) scenes.size() > 0)
        return scenes[currscene]->view;
    else
        return NULL;
}

Terrain * GLWidget::getTerrain()
{
    if((int) scenes.size() > 0)
        return scenes[currscene]->terrain;
    else
        return NULL;
}

TypeMap * GLWidget::getTypeMap(TypeMapType purpose)
{
    if((int) scenes.size() > 0)
    {
        return scenes[currscene]->maps[(int) purpose];
    }
    else
        return NULL;
}

PMrender::TRenderer * GLWidget::getRenderer()
{
    return renderer;
}

Skirts * GLWidget::getSkirts()
{
    if((int) scenes.size() > 0)
    {
        return scenes[currscene]->skirts;
    }
    else
        return NULL;
}

Plume * GLWidget::getPlume()
{
    if((int) scenes.size() > 0)
    {
        return scenes[currscene]->plume;
    }
    else
        return NULL;
}

AirLayers * GLWidget::getAirLayers()
{
    if((int) scenes.size() > 0)
    {
        return scenes[currscene]->airlayers;
    }
    else
        return NULL;
}

/*
BrushPalette * GLWidget::getPalette()
{
    return palette;
}*/

void GLWidget::refreshOverlay()
{
    renderer->updateTypeMapTexture(getTypeMap(scenes[currscene]->overlay), PMrender::TRenderer::typeMapInfo::PAINT, false);
    update();
}

void GLWidget::setOverlay(TypeMapType purpose)
{
    scenes[currscene]->overlay = purpose;
    renderer->updateTypeMapTexture(getTypeMap(scenes[currscene]->overlay), PMrender::TRenderer::typeMapInfo::PAINT, true);
    update();
}

TypeMapType GLWidget::getOverlay()
{
    return scenes[currscene]->overlay;
}

/*
float GLWidget::getRadius()
{
    return brushcursor.getRadius();
}

void GLWidget::setRadius(float rad)
{
    float brushradius = rad * getView()->getScaleConst();
    brushcursor.setRadius(brushradius);
    update();
}*/

void GLWidget::setMode(ControlMode mode)
{
    cmode = mode;
    
    switch(mode)
    {
        case ControlMode::VIEW:
               setOverlay(TypeMapType::DISPLAY);
               break;
        case ControlMode::PAINTBOIDS:
               setOverlay(TypeMapType::BOIDS);
               break;
        case ControlMode::PAINTTYPE:
               setOverlay(TypeMapType::CATEGORY);
               break;
        case ControlMode::PAINTSLOPE:
               setOverlay(TypeMapType::SLOPE);
               break;
        default:
               break;
    }
    update();
}

std::string GLWidget::get_dirprefix()
{
    // std::cout << "Datadir before fixing: " << datadir << std::endl;
    while (datadir.back() == '/')
        datadir.pop_back();

    // std::cout << "Datadir after fixing: " << datadir << std::endl;

    int slash_idx = datadir.find_last_of("/");
    std::string setname = datadir.substr(slash_idx + 1);
    std::string dirprefix = datadir + "/" + setname;
    return dirprefix;
}

void GLWidget::loadScene(int numframes, int start)
{
    std::cout << "Datadir before fixing: " << datadir << std::endl;
    while (datadir.back() == '/')
        datadir.pop_back();

    std::cout << "Datadir after fixing: " << datadir << std::endl;

    int slash_idx = datadir.find_last_of("/");
    std::string setname = datadir.substr(slash_idx + 1);
    std::string dirprefix = get_dirprefix();

    loadScene(dirprefix, numframes, start);
}

void GLWidget::loadScene(std::string dirprefix, int numframes, int start)
{
    std::string terfile = dirprefix+".elv";
    std::string skirtfile = datadir+"/";
    std::string plumefile = datadir+"/";
    std::string airlayerfile = datadir+"/";
    float tx, ty;
    int gx, gy;
    
    cerr << "LOAD SCENE" << endl;
    
    // load terrain
    currscene = 0;
    getTerrain()->loadElv(terfile, getTypeMap(TypeMapType::DISPLAY));
    getTerrain()->getTerrainDim(tx, ty);
    getTerrain()->getGridDim(gx, gy);
    cerr << "Elevation file loaded with dimensions " << tx << "m X " << ty << "m" << " and grid dimensions = " << gx << " X " << gy << endl;
    scf = getTerrain()->getMaxExtent();
    

    // getTypeMap(TypeMapType::DISPLAY)->matchDim(gx, gy);
    getTypeMap(TypeMapType::CATEGORY)->matchDim(gx, gy);
    getTypeMap(TypeMapType::BOIDS)->matchDim(gx, gy);
    
    getView()->setForcedFocus(getTerrain()->getFocus());
    getView()->setViewScale(getTerrain()->longEdgeDist());
    getView()->setDim(0.0f, 0.0f, (float) this->width(), (float) this->height());
    getTerrain()->calcMeanHeight();

    int fillval = 0;
    getTypeMap(TypeMapType::DISPLAY)->fill(fillval);
    // getTypeMap(TypeMapType::CATEGORY)->fill(fillval);
    // getTypeMap(TypeMapType::BOIDS)->fill(fillval);
    
   getSkirts()->setNumFrames(numframes);
   getSkirts()->loadSkirts(skirtfile, 4, 200);
    
    getPlume()->setNumFrames(numframes);
    getPlume()->loadPlume(plumefile, start);
    
//    getAirLayers()->setNumFrames(numframes);
//    getAirLayers()->loadAirLayers(airlayerfile, 5);
//    bindAirLayer();
  
    focuschange = true;
    setMode(ControlMode::VIEW);
    
    /*
    getPalette()->directSelect(getPalette()->firstSpeciesPalette());
    getPalette()->switchMode(ControlMode::PAINTTYPE);
    */
    
    getView()->apply();
    // setRadius(mainparent->getSliderValue());
    // mainparent->showPalette(true);
    
    update();

    std::cerr << "Scene Loading Finished" << std::endl;
}

void GLWidget::saveScene(std::string dirprefix)
{
    std::string terfile = dirprefix+".elv";

    // save terrain
    getTerrain()->saveElv(terfile);
}

bool GLWidget::readPaintMaps(std::string paintfile)
{
    int tp;
    int width, height;
    ifstream infile;

    infile.open((char *) paintfile.c_str(), ios_base::in);
    if(infile.is_open())
    {
        infile >> width >> height;
        getTypeMap(TypeMapType::CATEGORY)->matchDim(width, height);

        // write category map
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
            {
                infile >> tp;
                getTypeMap(TypeMapType::CATEGORY)->set(x,y,tp);
            }

        getTypeMap(TypeMapType::SLOPE)->matchDim(width, height);
        // write slope map
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
            {
                infile >> tp;
                getTypeMap(TypeMapType::SLOPE)->set(x,y,tp);
            }
        infile.close();
        return true;
    }
    else
    {
        cerr << "Error GLWidget::readPaintMaps: unable to open file" << paintfile << endl;
        return false;
    }
}

bool GLWidget::writePaintMaps(std::string paintfile)
{
    int width, height;
    ofstream outfile;

    outfile.open((char *) paintfile.c_str(), ios_base::out);
    if(outfile.is_open())
    {
        width = getTypeMap(TypeMapType::CATEGORY)->width();
        height = getTypeMap(TypeMapType::CATEGORY)->height();
        outfile << width << " " << height << endl;
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
            {
                outfile << getTypeMap(TypeMapType::CATEGORY)->get(x, y) << " ";
            }
        outfile << endl;
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
            {
                outfile << getTypeMap(TypeMapType::SLOPE)->get(x, y) << " ";
            }
        outfile << endl;
        outfile.close();
        return true;
    }
    else
    {
        cerr << "Error GLWidget::writePaintMaps:unable to open file " << paintfile << endl;
        return false;
    }
}

/*
void GLWidget::writePaintMap(std::string paintfile)
{
    getTypeMap(TypeMapType::CATEGORY)->saveToPaintImage(paintfile);
}
*/

void GLWidget::addScene()
{
    Scene * scene = new Scene();
    scene->view->setDim(0.0f, 0.0f, (float) this->width(), (float) this->height());

    scenes.push_back(scene);
    currscene = (int) scenes.size() - 1;
}

void GLWidget::setScene(int s)
{
    if(s >= 0 && s < (int) scenes.size())
    {
        currscene = s;
        getTerrain()->setBufferToDirty();
        refreshOverlay();
        update();
    }
}

void GLWidget::loadDecals()
{
    QImage decalImg, t;

    // load image
    
/*    if(!decalImg.load(QCoreApplication::applicationDirPath() + "/../../terviewer/Icons/manipDecals.png"))
        cerr << QCoreApplication::applicationDirPath().toUtf8().constData() << "/../../terviewer/Icons/manipDecals.png" << " not found" << endl;*/
    if(!decalImg.load(QCoreApplication::applicationDirPath() + "/../../test.ppm"))
            cerr << QCoreApplication::applicationDirPath().toUtf8().constData() << "/../../test.ppm" << " not found" << endl;

    // Qt prep image for OpenGL
    // QImage fixedImage(decalImg.width(), decalImg.height(), QImage::Format_ARGB32);
    QImage fixedImage(decalImg.width(), decalImg.height(), QImage::Format_RGB32);
    QPainter painter(&fixedImage);
    painter.setCompositionMode(QPainter::CompositionMode_Source);
    painter.fillRect(fixedImage.rect(), Qt::transparent);
    painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
    painter.drawImage( 0, 0, decalImg);
    painter.end();

    t = QGLWidget::convertToGLFormat( fixedImage );

    // renderer->bindDecals(t.width(), t.height(), t.bits());
    decalsbound = true;
}

int GLWidget::loadTypeMap(MapFloat * map, TypeMapType purpose)
{
    int numClusters = 0;

    switch(purpose)
    {
        case TypeMapType::BOIDS:
            break;
        case TypeMapType::CATEGORY:
            break;
        case TypeMapType::SLOPE:
            numClusters = getTypeMap(purpose)->convert(map, purpose, 1.0f);
            break;
        default:
            break;
    }
    return numClusters;
}

void GLWidget::positionPlantsToDisplay(int treespacing, int shrubspacing)
{
    int ax, ay, dx, dy;
    TypeMap * cat = getTypeMap(TypeMapType::CATEGORY);
    
    getTerrain()->getGridDim(dx, dy);
    
    trees.clear();
    shrubs.clear();
    
    // place trees
    for (int x = 0; x < dx; x+= treespacing)
        for (int y = 0; y < dy; y+= treespacing)
            if(cat->get(y, x) == (int) BrushType::FOREST)
            {
                Plant plnt;
                vpPoint p, gp;
                float h, r;
                glm::vec4 c;
                
                // create plant with random offsets to position, height, canopy radius and colour
                
                // jittered position on terrain with max separation according to spacing
                gp = getTerrain()->toWorld(x, y, 0.0f);
                vpPoint dpnt = vpPoint(gp.z, gp.y, gp.x);
                getTerrain()->drapePnt(dpnt, dpnt);
                gp.y = dpnt.y;
                float rndoffx = (nfield1->getNoise(gp) - 0.5f) * (float) treespacing * 0.9f;
                float rndoffy = (nfield2->getNoise(gp) - 0.5f) * (float) treespacing * 0.9f;
                ax = x + rndoffx;
                ay = y + rndoffy;
                
                // project onto terrain
                p = getTerrain()->toWorld(ax, ay, 0.0f);
                dpnt = vpPoint(p.z, p.y, p.x);
                getTerrain()->drapePnt(dpnt, dpnt);
                p.y = dpnt.y;
              
                float rndoff1 = nfield1->getNoise(gp)*0.3f-0.15f;
                c = glm::vec4(rndoff1, rndoff1, rndoff1, 1.0f); // randomly vary lightness of plant in [-0.15, 0.15]
                h = treeheight*(1.0f+rndoff1);
                float rndoff2 = nfield2->getNoise(gp)*0.3f-0.15f;
                r = treeradius*(0.85+rndoff2);
                
                plnt.pos = p;
                plnt.height = h;
                plnt.canopy = 2.0f * r;
                plnt.col = c;
                
                trees.push_back(plnt);
            }
    
    // place shrubs
    for (int x = 0; x < dx; x+= shrubspacing)
        for (int y = 0; y < dy; y+= shrubspacing)
            if(cat->get(y, x) == (int) BrushType::SHRUBS)
            {
                Plant plnt;
                vpPoint p, gp;
                float h, r;
                glm::vec4 c;
                
                // create plant with random offsets to position, height, canopy radius and colour
                
                // jittered position on terrain with max separation according to spacing
                gp = getTerrain()->toWorld(x, y, 0.0f);
                vpPoint dpnt = vpPoint(gp.z, gp.y, gp.x);
                getTerrain()->drapePnt(dpnt, dpnt);
                gp.y = dpnt.y;
                float rndoffx = (nfield2->getNoise(gp) - 0.5f) * (float) shrubspacing * 1.1f;
                float rndoffy = (nfield1->getNoise(gp) - 0.5f) * (float) shrubspacing * 1.1f;
                ax = x + rndoffx;
                ay = y + rndoffy;
                
                // project onto terrain
                p = getTerrain()->toWorld(ax, ay, 0.0f);
                dpnt = vpPoint(p.z, p.y, p.x);
                getTerrain()->drapePnt(dpnt, dpnt);
                p.y = dpnt.y;
              
                float rndoff1 = nfield2->getNoise(gp)*0.2f-0.1f;
                c = glm::vec4(rndoff1, rndoff1, rndoff1, 1.0f); // randomly vary lightness of plant in [-0.15, 0.15]
                h = shrubheight*(1.0f+rndoff1);
                float rndoff2 = nfield1->getNoise(gp)*0.3f-0.15f;
                r = shrubradius*(0.85+rndoff2);
                
                plnt.pos = p;
                plnt.height = h;
                plnt.canopy = 2.0f * r;
                plnt.col = c;
                
                shrubs.push_back(plnt);
            }
}

void GLWidget::bindAirLayer()
{
    QImage & timg = getAirLayers()->getCurrentLayerImage();
    renderer->bindDecals(timg.width(), timg.height(), timg.bits());
    decalsbound = true;
    update();
}

void GLWidget::genAirLayerShape()
{
    glm::mat4 idt, tfm;
    glm::vec3 rotx;
 
    // plane
    rotx = glm::vec3(1.0f, 0.0f, 0.0f);
    idt = glm::mat4(1.0f);
    tfm = glm::rotate(idt, glm::radians(-90.0f), rotx);
    airlayershape.genPlane(tfm);
}

void GLWidget::genSkirtShape()
{
    glm::mat4 idt;
 
    // sphere
    idt = glm::mat4(1.0f);
#ifdef HIGHRES
    skirtshape.genSphere(1.0f, 20, 20, idt);
#endif
#ifdef LOWRES
    skirtshape.genSphere(1.0f, 5, 5, idt);
#endif
}

void GLWidget::genPlumeShape()
{
    glm::mat4 idt;
 
    // sphere
    idt = glm::mat4(1.0f);
#ifdef HIGHRES
    plumeshape.genSphere(1.0f, 20, 20, idt);
#endif
#ifdef LOWRES
    plumeshape.genSphere(1.0f, 5, 5, idt);
#endif
}

void GLWidget::genTreeShape(float trunkheight, float trunkradius)
{
    glm::mat4 idt, tfm;
    glm::vec3 trs, rotx;
    float canopyheight;

    rotx = glm::vec3(1.0f, 0.0f, 0.0f);
    canopyheight = 1.0f - trunkheight;

    // trunk - cylinder
    idt = glm::mat4(1.0f);
    tfm = glm::rotate(idt, glm::radians(-90.0f), rotx);
    // extra 0.1 to trunk is to ensure that trunk joins tree properly
    treeshape.genCappedCylinder(trunkradius, trunkradius, trunkheight+0.1f, 3, 1, tfm, false);

    // canopy - cone
    idt = glm::mat4(1.0f);
    trs = glm::vec3(0.0f, trunkheight, 0.0f);
    tfm = glm::translate(idt, trs);
    tfm = glm::rotate(tfm, glm::radians(-90.0f), rotx);
#ifdef HIGHRES
    treeshape.genCappedCone(0.5f, canopyheight, 20, 1, tfm, false);
#endif
#ifdef LOWRES
    treeshape.genCappedCone(0.5f, canopyheight, 6, 1, tfm, false);
#endif
}

void GLWidget::genShrubShape(float trunkheight, float trunkradius)
{
    glm::mat4 idt, tfm;
    glm::vec3 trs, rotx;
    float canopyheight;

    rotx = glm::vec3(1.0f, 0.0f, 0.0f);
    canopyheight = 1.0f - trunkheight;

    // trunk - cylinder
    idt = glm::mat4(1.0f);
    tfm = glm::rotate(idt, glm::radians(-90.0f), rotx);
    // extra 0.1 to trunk is to ensure that trunk joins tree properly
    shrubshape.genCappedCylinder(trunkradius, trunkradius, trunkheight+0.1f, 3, 1, tfm, false);

    // canopy - sphere
    idt = glm::mat4(1.0f);
    trs = glm::vec3(0.0f, trunkheight+canopyheight/2.0f, 0.0f);
    tfm = glm::translate(idt, trs);
    tfm = glm::scale(tfm, glm::vec3(1.0, canopyheight, 1.0f)); // make sure tree fills 1.0f on a side bounding box
    tfm = glm::rotate(tfm, glm::radians(-90.0f), rotx);

#ifdef HIGHRES
    shrubshape.genSphere(0.5f, 20, 20, tfm);
#endif
#ifdef LOWRES
    shrubshape.genSphere(0.5f, 6, 6, tfm);
#endif
    
}

void GLWidget::drawAirLayers(std::vector<ShapeDrawData> &drawParams)
{
    ShapeDrawData sdd;
    
    std::vector<glm::mat4> xforms; // transformation to be applied to each instance
    std::vector<glm::vec4> colvars; // colour variation to be applied to each instance
    
    if(getAirLayers()->getNumFrames() > 0)
    {
        std::vector<AirLayer> &airlayers = getAirLayers()->getCurrentPlanes();
        
        xforms.clear();
        colvars.clear();
        
        // Gather plume instances
        for(auto &lyr: airlayers)
        {
            if(lyr.active)
            {
                // setup transformation for individual plant, including scaling and translation
                glm::mat4 idt, tfm;
                glm::vec3 trs, sc;
                
                idt = glm::mat4(1.0f);
                trs = glm::vec3(30000.0f, lyr.altitude, 30000.0f);
                tfm = glm::translate(idt, trs);
                sc = glm::vec3(60000.0f, 100.0f, 60000.0f); // TO DO JG - get terrain extent
                tfm = glm::scale(tfm, sc);
                xforms.push_back(tfm); // affine transform for placement
                
                // glm::vec4 ashcol = {1.0f - sph.ashloss, 1.0f - sph.ashloss, 1.0f - sph.ashloss, 1.0f};
                glm::vec4 ashcol = {0.0f, 0.0f, 0.0f, 1.0f};
                colvars.push_back(ashcol); // colour variation due to ash loss
            }
        }
    }
    
    airlayershape.removeAllInstances();
    airlayershape.bindInstances(&xforms, &colvars);
    
    // prep drawing
    sdd = airlayershape.getDrawParameters();
    sdd.current = false;
    drawParams.push_back(sdd);
}

void GLWidget::drawSkirts(std::vector<ShapeDrawData> &drawParams)
{
    ShapeDrawData sdd;
    
    std::vector<glm::mat4> xforms; // transformation to be applied to each instance
    std::vector<glm::vec4> colvars; // colour variation to be applied to each instance
    
    // cerr << "draw skirts" << endl;
    // cerr << "frames = " << getSkirts()->getNumFrames() << endl;
    // cerr << "layers = " << getSkirts()->getNumLayers() << endl;
    
    if(getSkirts()->getNumFrames() > 0)
    {
        xforms.clear();
        colvars.clear();
        
        for(int lyr = 0; lyr < getSkirts()->getNumLayers(); lyr++)
        {
            for(int s = 0; s < getSkirts()->getNumSamples(lyr); s++)
            {
                // setup transformation for skirt, including scaling and translation
                glm::mat4 idt, tfm;
                glm::vec3 trs, sc;
                vpPoint pos;
                float val, rad = 50.0f;
                
                // cerr << "lyr = " << lyr << " sample = " << s << " of " << getSkirts()->getNumSamples(lyr) << endl;
                getSkirts()->getFrameSample(lyr, s, pos, val);
                
                if(val > 0.01) // only rendered condensed samples
                {
                    // cerr << sph.pos.x << " " << sph.pos.y << " " << sph.pos.z << endl;
                    idt = glm::mat4(1.0f);
                    trs = glm::vec3(pos.x, pos.y, pos.z);
                    tfm = glm::translate(idt, trs);
                    sc = glm::vec3(rad, rad, rad);
                    tfm = glm::scale(tfm, sc);
                    xforms.push_back(tfm); // affine transform for placement
                    
                    glm::vec4 chgcol;
                    /*
                    if(val > 5.0f)
                        chgcol = {0.5f, 0.0f, -0.5f, 1.0f};
                    else
                        chgcol = {0.0f, 0.0f, 0.0f, 1.0f};
                    */
                    switch(lyr)
                    {
                        case 0:
                            chgcol = {-1.0f, -1.0f, 0.0f, 1.0f};
                            break;
                        case 1:
                            chgcol = {0.0f, -1.0f, -1.0f, 1.0f};
                            break;
                        case 2:
                            chgcol = {-1.0f, 1.0f, -1.0f, 1.0f};
                            break;
                        default:
                            break;
                    }
                    colvars.push_back(chgcol); // colour variation due to condensation
                }
                // cerr << "x = " << x << " y = " << y << endl;
            }
        }
    }
        
    skirtshape.removeAllInstances();
    skirtshape.bindInstances(&xforms, &colvars);
    
    // prep drawing
    sdd = skirtshape.getDrawParameters();
    sdd.current = false;
    drawParams.push_back(sdd);
}

void GLWidget::drawPlume(std::vector<ShapeDrawData> &drawParams)
{
    ShapeDrawData sdd;
    
    std::vector<glm::mat4> xforms; // transformation to be applied to each instance
    std::vector<glm::vec4> colvars; // colour variation to be applied to each instance
    
    if(getPlume()->getNumFrames() > 0)
    {
        std::vector<PlumeSphere> &plumespheres = getPlume()->getCurrentPlume();
        
        xforms.clear();
        colvars.clear();
        
        // cerr << "Current plume size " << (int) plumespheres.size() << endl;;
        // Gather plume instances
        for(auto &sph: plumespheres)
        {
            // setup transformation for individual plant, including scaling and translation
            glm::mat4 idt, tfm;
            glm::vec3 trs, sc;
            
            // cerr << sph.pos.x << " " << sph.pos.y << " " << sph.pos.z << endl;
            idt = glm::mat4(1.0f);
            trs = glm::vec3(sph.pos.x, sph.pos.y, sph.pos.z);
            tfm = glm::translate(idt, trs);
            sc = glm::vec3(sph.radius, sph.radius, sph.radius);
            tfm = glm::scale(tfm, sc);
            xforms.push_back(tfm); // affine transform for placement
            
            if(plumecharge) // colour according to charge
            {
                // high charge closer to red
                /*
                float c = fmaxf(0.0, sph.charge - 100.0f) / 100.0f;
                glm::vec4 chgcol = {1.0f, 1.0f-c, 1.0f-c, 0.0f};
                 */
                glm::vec4 chgcol;
                switch(sph.sphtype)
                {
                    case 0:
                        chgcol = {1.0f, 0.0f, 0.0f, 1.0f};
                        break;
                    case 1:
                        chgcol = {0.0f, 1.0f, 0.0f, 1.0f};
                        break;
                    case 2:
                        chgcol = {0.0f, 0.0f, 1.0f, 1.0f};
                        break;
                    case 3:
                        chgcol = {1.0f, 0.0f, 1.0f, 1.0f};
                        break;
                }
                colvars.push_back(chgcol); // colour variation due to charge
                
            }
            else // colour according to density
            {
                // high density closer to black
                float d = fmaxf(0.0, 1.5f - sph.density) / 1.75f;
                // float d = 0.0f;
                glm::vec4 ashcol = {d, d, d, 0.0f};
                colvars.push_back(ashcol); // colour variation due to ash loss
            }
          
        }
    }
    
    plumeshape.removeAllInstances();
    plumeshape.bindInstances(&xforms, &colvars);
    
    // prep drawing
    sdd = plumeshape.getDrawParameters();
    sdd.current = false;
    drawParams.push_back(sdd);
}

void GLWidget::drawPlants(std::vector<ShapeDrawData> &drawParams)
{
    ShapeDrawData sdd;
    
    std::vector<glm::mat4> xforms; // transformation to be applied to each instance
    std::vector<glm::vec4> colvars; // colour variation to be applied to each instance
    
    xforms.clear();
    colvars.clear();
    
    // Gather tree instances
    for(auto &plnt: trees)
    {
        // setup transformation for individual plant, including scaling and translation
        glm::mat4 idt, tfm;
        glm::vec3 trs, sc;
        // glm::vec3 rotate_axis = glm::vec3(0.0f, 1.0f, 0.0f);
        vpPoint loc = plnt.pos;
        // GLfloat rotate_rad;
        
        idt = glm::mat4(1.0f);
        trs = glm::vec3(loc.x, loc.y, loc.z);
        tfm = glm::translate(idt, trs);
        sc = glm::vec3(plnt.canopy, plnt.height, plnt.canopy);
        tfm = glm::scale(tfm, sc);
        // tfm = glm::rotate(tfm, rotate_rad, rotate_axis);
        xforms.push_back(tfm); // affine transform for placement
        colvars.push_back(plnt.col); // colour variation
    }
    
    treeshape.removeAllInstances();
    treeshape.bindInstances(&xforms, &colvars);
    
    // prep drawing
    sdd = treeshape.getDrawParameters();
    sdd.current = false;
    drawParams.push_back(sdd);
    
    xforms.clear();
    colvars.clear();
    
    // Gather shrub instances
    for(auto &plnt: shrubs)
    {
        // setup transformation for individual plant, including scaling and translation
        glm::mat4 idt, tfm;
        glm::vec3 trs, sc;
        // glm::vec3 rotate_axis = glm::vec3(0.0f, 1.0f, 0.0f);
        vpPoint loc = plnt.pos;
        // GLfloat rotate_rad;
        
        idt = glm::mat4(1.0f);
        trs = glm::vec3(loc.x, loc.y, loc.z);
        tfm = glm::translate(idt, trs);
        sc = glm::vec3(plnt.canopy, plnt.height, plnt.canopy);
        tfm = glm::scale(tfm, sc);
        // tfm = glm::rotate(tfm, rotate_rad, rotate_axis);
        xforms.push_back(tfm); // affine transform for placement
        colvars.push_back(plnt.col); // colour variation
    }
    
    shrubshape.removeAllInstances();
    shrubshape.bindInstances(&xforms, &colvars);
    
    // prep drawing
    sdd = shrubshape.getDrawParameters();
    sdd.current = false;
    drawParams.push_back(sdd);
}

void GLWidget::initializeGL()
{
    // get context opengl-version
    qDebug() << "Widget OpenGl: " << format().majorVersion() << "." << format().minorVersion();
    qDebug() << "Context valid: " << context()->isValid();
    qDebug() << "Really used OpenGl: " << context()->format().majorVersion() << "." <<
              context()->format().minorVersion();
    qDebug() << "OpenGl information: VENDOR:       " << (const char*)glGetString(GL_VENDOR);
    qDebug() << "                    RENDERDER:    " << (const char*)glGetString(GL_RENDERER);
    qDebug() << "                    VERSION:      " << (const char*)glGetString(GL_VERSION);
    qDebug() << "                    GLSL VERSION: " << (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);

    QGLFormat glFormat = QGLWidget::format();
    if ( !glFormat.sampleBuffers() )
        qWarning() << "Could not enable sample buffers";

    qglClearColor(qtWhite.light());

    int mu;
    glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &mu);
    cerr << "max texture units = " << mu << endl;

    // *** PM Render code - start ***

    // To use basic shading: PMrender::TRenderer::BASIC
    // To use radianvce scaling: PMrender::TRenderer::RADIANCE_SCALING

    PMrender::TRenderer::terrainShadingModel sMod = PMrender::TRenderer::RADIANCE_SCALING;

    // set terrain shading model
    renderer->setTerrShadeModel(sMod);

    // set up light
    Vector dl = Vector(0.6f, 1.0f, 0.6f);
    dl.normalize();

    GLfloat pointLight[3] = { 0.5, 5.0, 7.0}; // side panel + BASIC lighting
    GLfloat dirLight0[3] = { dl.i, dl.j, dl.k}; // for radiance lighting
    GLfloat dirLight1[3] = { -dl.i, dl.j, -dl.k}; // for radiance lighting

    renderer->setPointLight(pointLight[0],pointLight[1],pointLight[2]);
    renderer->setDirectionalLight(0, dirLight0[0], dirLight0[1], dirLight0[2]);
    renderer->setDirectionalLight(1, dirLight1[0], dirLight1[1], dirLight1[2]);

    // initialise renderer/compile shaders
    renderer->initShaders();

    // set other render parameters
    // can set terrain colour for radiance scaling etc - check trenderer.h

    // terrain contours
    renderer->drawContours(false);
    renderer->drawGridlines(false);

    // turn on terrain type overlay (off by default); NB you can stil call methods to update terrain type,
    renderer->useTerrainTypeTexture(true);
    renderer->useConstraintTypeTexture(false);

    // use manipulator textures (decal'd)
    renderer->textureManipulators(true);

    // *** PM Render code - end ***

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_DEPTH_CLAMP);
    glEnable(GL_TEXTURE_2D);

    loadDecals();
    paintGL();
}

void GLWidget::paintGL()
{
    vpPoint mo;
    glm::mat4 tfm, idt;
    glm::vec3 trs, rot;
    uts::vector<ShapeDrawData> drawParams; // to be passed to terrain renderer
    Shape shape;  // geometry for focus indicator
    std::vector<glm::mat4> sinst;
    std::vector<glm::vec4> cinst;

    Timer t;

    if(active)
    {
        t.start();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // note: bindinstances will not work on the first render pass because setup time is needed

        if(focuschange && focusviz)
        {
            ShapeDrawData sdd;
            float scale;

            GLfloat manipCol[] = {0.325f, 0.235f, 1.0f, 1.0f};

            // create shape
            shape.clear();
            shape.setColour(manipCol);

            // place vertical cylinder at view focus
            mo = getView()->getFocus();
            scale = getView()->getScaleFactor();
            idt = glm::mat4(1.0f);
            trs = glm::vec3(mo.x, mo.y, mo.z);
            rot = glm::vec3(1.0f, 0.0f, 0.0f);
            tfm = glm::translate(idt, trs);
            tfm = glm::rotate(tfm, glm::radians(-90.0f), rot);
            shape.genCappedCylinder(scale*armradius, 1.5f*scale*armradius, scale*(manipheight-manipradius), 40, 10, tfm, false);
            if(shape.bindInstances(&sinst, &cinst)) // passing in an empty instance will lead to one being created at the origin
            {
                sdd = shape.getDrawParameters();
                sdd.current = false;
                drawParams.push_back(sdd);
            }
        }

        if(cmode == ControlMode::PLACEHERD || cmode == ControlMode::PAINTBOIDS || cmode == ControlMode::PAINTTYPE || cmode == ControlMode::PAINTSLOPE)
        {
            // ShapeDrawData sdd;

            // cerr << "Set Brush" << endl;
            // cerr << "draw type = " << (int) palette->getDrawType() << endl;
            
            /*
            brushcursor.setBrushColour(getTypeMap(TypeMapType::CATEGORY)->getColour((int) palette->getDrawType()));
            brushcursor.shape.clear();
            brushcursor.genBrushRing(getView(), getTerrain(), getRadius(), false);
            if(brushcursor.shape.bindInstances(&sinst, &cinst))
            {
                sdd = brushcursor.shape.getDrawParameters();
                sdd.current = false;
                drawParams.push_back(sdd);
            }*/

            /*
            brushcursor.debugshape.clear();
            if(brushcursor.debugshape.bindInstances(&sinst, &cinst))
            {
                sdd = brushcursor.debugshape.getDrawParameters();
                sdd.current = false;
                drawParams.push_back(sdd);
            }
            */
            // cerr << "Ring radius = " << getRadius() << endl;
        }

        // drawAirLayers(drawParams);
        if(getPlumeVisible())
            drawPlume(drawParams);
        if(getSkirtsVisible())
            drawSkirts(drawParams);

        // pass in draw params for objects
        renderer->setConstraintDrawParams(drawParams);

        // draw terrain and plants
        getTerrain()->updateBuffers(renderer); 

        if(focuschange)
            renderer->updateTypeMapTexture(getTypeMap(getOverlay())); // only necessary if the texture is changing dynamically
        renderer->draw(getView());

        t.stop();

        if(timeron)
            cerr << "rendering = " << t.peek() << " fps = " << 1.0f / t.peek() << endl;
    }
}

void GLWidget::resizeGL(int width, int height)
{
    // TO DO: fix resizing
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, width, height);

    cerr << "RESIZE " << this->width() << " X " << this->height() << endl;
    // apply to all views
    for(int i = 0; i < (int) scenes.size(); i++)
    {
        scenes[i]->view->setDim(0.0f, 0.0f, (float) this->width(), (float) this->height());
        scenes[i]->view->apply();
    }
}


void GLWidget::keyPressEvent(QKeyEvent *event)
{
    /*
    if(event->key() == Qt::Key_A) // 'A' for animated spin around center point of terrain
    {
        getView()->startSpin();
        rtimer->start(20);
    }
    */
    if(event->key() == Qt::Key_A) // 'A' toggle visibility of planes
    {
        getAirLayers()->setActiveStatus(!getAirLayers()->getActiveStatus());
        update();
    }
    /*
    if(event->key() == Qt::Key_C) // 'C' to capture an image
    {
        QImage * capImg = new QImage();
        // QImage fixedImage(decalImg.width(), decalImg.height(), QImage::Format_ARGB32);
        unscaledScreenCapture(capImg);
        
        std::string filename = datadir+"cap.png";
        capImg->save(QString::fromStdString(filename));
        cerr << "Screen Capture Taken" << endl;
    }*/
    if(event->key() == Qt::Key_B) // 'B' for condensed water layer view
    {
        getAirLayers()->setMode(VolcanoMode::CONDENSED);
        bindAirLayer();
    }
    if(event->key() == Qt::Key_C) // 'C' for convection layer view
    {
        getAirLayers()->setMode(VolcanoMode::CONVECTION);
        bindAirLayer();
    }
    /*
    if(event->key() == Qt::Key_F) // 'F' to toggle focus stick visibility
    {
        if(focusviz)
            focusviz = false;
        else
            focusviz = true;
        update();
    }*/
    if(event->key() == Qt::Key_F) // 'F' for ash layer view
    {
        getAirLayers()->setMode(VolcanoMode::ASH);
        bindAirLayer();
    }
    if(event->key() == Qt::Key_G) // 'G' toggle between charge and density display
    {
        setPlumeCharge(!getPlumeCharge());
        update();
    }
    if(event->key() == Qt::Key_I) // 'I' move up one weather layer
    {
        int lyr = getAirLayers()->getLayer();
        if(lyr < getAirLayers()->getNumLayers()-1)
        {
            getAirLayers()->getCurrentPlanes()[lyr].active = false;
            getAirLayers()->getCurrentPlanes()[lyr+1].active = true;
            getAirLayers()->setLayer(getAirLayers()->getLayer()+1);
            bindAirLayer();
        }
        cerr << "Layer = " << getAirLayers()->getLayer() << endl;
    }
    if(event->key() == Qt::Key_M) // 'M' for moisture layer view
    {
        getAirLayers()->setMode(VolcanoMode::MOISTURE);
        bindAirLayer();
    }
    if(event->key() == Qt::Key_O) // 'O' step back volcano frame
    {
        getPlume()->setFrame(getPlume()->getFrame()-1);
        getSkirts()->setFrame(getSkirts()->getFrame()-1);
        // getAirLayers()->setFrame(getAirLayers()->getFrame()-1);
        // bindAirLayer();
        update();
        cerr << "Sim Step = " << getPlume()->getFrame() << endl;
    }
    if(event->key() == Qt::Key_P) // 'P' step forward volcano frame
    {
        getPlume()->setFrame(getPlume()->getFrame()+1);
        getSkirts()->setFrame(getSkirts()->getFrame()+1);
        // getAirLayers()->setFrame(getAirLayers()->getFrame()+1);
        // bindAirLayer();
        update();
        cerr << "Sim Step = " << getPlume()->getFrame() << endl;
    }
    if(event->key() == Qt::Key_R) // 'R' to restore saved view state
    {
        cerr << "R pressed" << endl;
        std::string vfile = datadir+"/view.txt";
        cerr << vfile << endl;
        getView()->load(vfile.c_str());
        getView()->apply();
        update();
    }
    /*
    if(event->key() == Qt::Key_S) // 'S' to save view state
    {
        cerr << "S pressed" << endl;
        std::string vfile = datadir+"/view.txt";
        cerr << vfile << endl;
        getView()->save(vfile.c_str());
        getView()->apply();
        update();
    }*/
    if(event->key() == Qt::Key_S) // 'S' toggle visibility of skirts
    {
        setSkirtsVisible(!getSkirtsVisible());
        update();
    }
  
    /*
    if(event->key() == Qt::Key_T) // 'T' for top-down view
    {
        getTerrain()->setMidFocus();
        getView()->setForcedFocus(getTerrain()->getFocus());
        getView()->topdown();
        update();
    }
     */
    if(event->key() == Qt::Key_T) // 'T' for temperature layer view
    {
        getAirLayers()->setMode(VolcanoMode::TEMPERATURE);
        bindAirLayer();
    }
    if(event->key() == Qt::Key_U) // 'U' move down one weather layer
    {
        int lyr = getAirLayers()->getLayer();
        if(lyr > 0)
        {
            getAirLayers()->getCurrentPlanes()[lyr].active = false;
            getAirLayers()->getCurrentPlanes()[lyr-1].active = true;
            getAirLayers()->setLayer(getAirLayers()->getLayer()-1);
            bindAirLayer();
        }
        cerr << "Layer = " << getAirLayers()->getLayer() << endl;
    }
    if(event->key() == Qt::Key_V) // 'V' toggle visibility of plume
    {
        setPlumeVisible(!getPlumeVisible());
        update();
    }
    if(event->key() == Qt::Key_X) // 'X' for xvelocity layer view
    {
        getAirLayers()->setMode(VolcanoMode::VELOCITYX);
        bindAirLayer();
    }
    if(event->key() == Qt::Key_Y) // 'Y' for yvelocity layer view
    {
        getAirLayers()->setMode(VolcanoMode::VELOCITYY);
        bindAirLayer();
    }
  
    /*
    if(event->key() == Qt::Key_V) // 'V' toggle to viewing mode
    {
        setMode(ControlMode::VIEW);
        getPalette()->switchMode(ControlMode::VIEW);
        getView()->apply();
        mainparent->showPalette(false);
        // qApp->setOverrideCursor( QCursor( Qt::ArrowCursor ) );
        update();
    }*/
    
    /*
    if(event->key() == Qt::Key_1) // activate/deactivate layer 1
    {
        if(getAirLayers()->getNumLayers() > 0)
        {
            getAirLayers()->getCurrentPlanes()[0].active = !getAirLayers()->getCurrentPlanes()[0].active;
            update();
        }
    }
    if(event->key() == Qt::Key_2) // activate/deactivate layer 2
    {
        if(getAirLayers()->getNumLayers() > 1)
        {
            getAirLayers()->getCurrentPlanes()[1].active = !getAirLayers()->getCurrentPlanes()[1].active;
            update();
        }
    }
    if(event->key() == Qt::Key_3) // activate/deactivate layer 3
    {
        if(getAirLayers()->getNumLayers() > 2)
        {
            getAirLayers()->getCurrentPlanes()[2].active = !getAirLayers()->getCurrentPlanes()[2].active;
            update();
        }
    }
    if(event->key() == Qt::Key_4) // activate/deactivate layer 4
    {
        if(getAirLayers()->getNumLayers() > 3)
        {
            getAirLayers()->getCurrentPlanes()[3].active = !getAirLayers()->getCurrentPlanes()[3].active;
            update();
        }
    }
    if(event->key() == Qt::Key_5) // activate/deactivate layer 5
    {
        if(getAirLayers()->getNumLayers() > 4)
        {
            getAirLayers()->getCurrentPlanes()[4].active = !getAirLayers()->getCurrentPlanes()[4].active;
            update();
        }
    }
    if(event->key() == Qt::Key_6) // activate/deactivate layer 6
    {
        if(getAirLayers()->getNumLayers() > 5)
        {
            getAirLayers()->getCurrentPlanes()[5].active = !getAirLayers()->getCurrentPlanes()[5].active;
            update();
        }
    }
    if(event->key() == Qt::Key_7) // activate/deactivate layer 7
    {
        if(getAirLayers()->getNumLayers() > 6)
        {
            getAirLayers()->getCurrentPlanes()[6].active = !getAirLayers()->getCurrentPlanes()[6].active;
            update();
        }
    }
    if(event->key() == Qt::Key_8) // activate/deactivate layer 8
    {
        if(getAirLayers()->getNumLayers() > 7)
        {
            getAirLayers()->getCurrentPlanes()[7].active = !getAirLayers()->getCurrentPlanes()[7].active;
            update();
        }
    }
    if(event->key() == Qt::Key_9) // activate/deactivate layer 9
    {
        if(getAirLayers()->getNumLayers() > 8)
        {
            getAirLayers()->getCurrentPlanes()[8].active = !getAirLayers()->getCurrentPlanes()[8].active;
            update();
        }
    }
    if(event->key() == Qt::Key_0) // activate/deactivate layer 10
    {
        if(getAirLayers()->getNumLayers() > 9)
        {
            getAirLayers()->getCurrentPlanes()[9].active = !getAirLayers()->getCurrentPlanes()[9].active;
            update();
        }
    }
     */
    
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    float nx, ny;
    vpPoint pnt;
    
    int x = event->x(); int y = event->y();
    float W = (float) width(); float H = (float) height();

    update(); // ensure this viewport is current for unproject

    // control view orientation with right mouse button or ctrl/alt modifier key and left mouse
    if(!viewlock && (event->modifiers() == Qt::MetaModifier || event->modifiers() == Qt::AltModifier || event->buttons() == Qt::RightButton))
    {
        // arc rotate in perspective mode
  
        // convert to [0,1] X [0,1] domain
        nx = (2.0f * (float) x - W) / W;
        ny = (H - 2.0f * (float) y) / H;
        lastPos = event->pos();
        getView()->startArcRotate(nx, ny);
        viewing = true;
    }
    lastPos = event->pos();
}

void GLWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    // set the focus for arcball rotation
    // pick point on terrain or zero plane if outside the terrain bounds
    vpPoint pnt;
    int sx, sy;
    
    sx = event->x(); sy = event->y();
    if(!viewlock && ((event->modifiers() == Qt::MetaModifier && event->buttons() == Qt::LeftButton) || (event->modifiers() == Qt::AltModifier && event->buttons() == Qt::LeftButton) || event->buttons() == Qt::RightButton))
    {
        getView()->apply();

        /*
        int port[4];
        getRenderer()->getViewPortDim(port);
        getView()->projectingPointPort(sx, sy, pnt, port[0], port[1], port[2], port[3]);*/

        // trackpnt = pnt;
        // trackpnt = view->getCOP();

        if(getTerrain()->pick(sx, sy, getView(), pnt))
        {
            if(!decalsbound)
                loadDecals();
            vpPoint pickpnt = pnt;
            getView()->setAnimFocus(pickpnt);
            getTerrain()->setFocus(pickpnt);
            // cerr << "Pick Point = " << pickpnt.x << ", " << pickpnt.y << ", " << pickpnt.z << endl;
            focuschange = true; focusviz = true;
            atimer->start(10);
        }
        // ignores pick if terrain not intersected, should possibly provide error message to user
    }
}

void GLWidget::pickInfo(int x, int y)
{
   std::string catName;

   cerr << endl;
   cerr << "*** PICK INFO ***" << endl;
   cerr << "location: " << x << ", " << y << endl;
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    viewing = false;

    /*
    if(event->button() == Qt::LeftButton && cmode == ControlMode::VIEW) // info on terrain cell
    {
        vpPoint pnt;
        int sx, sy;

        sx = event->x(); sy = event->y();

        if(getTerrain()->pick(sx, sy, getView(), pnt))
        {
            int x, y;
            getTerrain()->toGrid(pnt, x, y);
            pickInfo(x, y);
        }
    }
     */
   
    /*
    if(event->button() == Qt::LeftButton)
    {
        brush.finStroke();
        // getTypeMap(TypeMapType::CATEGORY)->setRegion(getTerrain()->coverRegion());
        renderer->updateTypeMapTexture(getTypeMap(getOverlay()));
    }*/
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    float nx, ny, W, H;
    // int port[4];

    int x = event->x();
    int y = event->y();

    W = (float) width();
    H = (float) height();

    // control view orientation with right mouse button or ctrl modifier key and left mouse
    if(!viewlock && ((event->modifiers() == Qt::MetaModifier && event->buttons() == Qt::LeftButton) || (event->modifiers() == Qt::AltModifier && event->buttons() == Qt::LeftButton) || event->buttons() == Qt::RightButton))
    {
        // convert to [0,1] X [0,1] domain
        nx = (2.0f * (float) x - W) / W;
        ny = (H - 2.0f * (float) y) / H;
        getView()->arcRotate(nx, ny);
        update();
        lastPos = event->pos();
    }

    //  if(!(event->buttons() == Qt::AllButtons) && cmode == ControlMode::PAINTTYPE) // show brush outline, whether or not the mouse is down
    if(cmode == ControlMode::PLACEHERD || cmode == ControlMode::PAINTBOIDS || cmode == ControlMode::PAINTTYPE || cmode == ControlMode::PAINTSLOPE)
    {
        // show brush radius
        // renderer->getViewPortDim(port);
        // brushcursor.cursorUpdate(getView(), getTerrain(), x, y, port);
        update();
    }
}

void GLWidget::wheelEvent(QWheelEvent * wheel)
{
    float del;
 
    QPoint pix = wheel->pixelDelta();
    QPoint deg = wheel->angleDelta();

    if(!viewlock)
    {
        if(!pix.isNull()) // screen resolution tracking, e.g., from magic mouse
        {
            del = (float) pix.y() * 10.0f;
            getView()->incrZoom(del);
            update();

        }
        else if(!deg.isNull()) // mouse wheel instead
        {
            del = (float) deg.y() * 2.5f;
            getView()->incrZoom(del);
            update();
        }
    }
}

void GLWidget::animUpdate()
{
    if(getView()->animate())
        update();
}

void GLWidget::rotateUpdate()
{
    if(getView()->spin())
        update();
}
