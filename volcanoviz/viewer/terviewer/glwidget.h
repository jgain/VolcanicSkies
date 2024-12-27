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

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include "glheaders.h" // Must be included before QT opengl headers
#include <QGLWidget>
#include <QLabel>
#include <QTimer>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QPushButton>
#include <list>
#include "common/basic_types.h"

#include "view.h"
#include "typemap.h"
#include "stroke.h"
#include "shape.h"
#include "dice_roller.h"
// #include "common/basic_types.h"
#include "kdtree.h"
#include "palette.h"
#include "noisefield.h"
#include "volcano.h"

/*
#include <QWidget>
#include <QImage>
#include <common/map.h>
*/

#define LOWRES


//! [0]

const float initmint = 0.0f;
const float initmaxt = 40.0f;
const float mtoft = 3.28084f;

const float treeheight = 5.0f;
const float shrubheight = 0.8f;
const float treeradius = 1.5f;
const float shrubradius = 0.6f;

class Window;

class Scene
{
public:

    View * view;
    Terrain * terrain;
    Plume * plume;
    Skirts * skirts;
    AirLayers * airlayers;
    TypeMap * maps[(int) TypeMapType::TMTEND];
    TypeMapType overlay;    //< currently active overlay texture

    Scene();

    ~Scene();
};

struct Plant
{
    vpPoint pos;    //< position on terrain in world units
    float height;   //< plant height in metres
    float canopy;   //< canopy radius in metres
    glm::vec4 col;  //< colour variation randomly assigned to plant
};

class Window;

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:

    GLWidget(const QGLFormat& format, string datadir, QWidget *parent = 0);
    
    ~GLWidget();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    /**
     * capture the framebuffer as an image
     * @param capImg    framebuffer is written to this image
     * @param capSize   the image is scaled by linear interpolation to this size
     */
    void screenCapture(QImage * capImg, QSize capSize);
    void unscaledScreenCapture(QImage * capImg);

    /// getters for currently active view, terrain, typemaps, renderer, ecosystem
    View * getView();
    Terrain * getTerrain();
    TypeMap * getTypeMap(TypeMapType purpose);
    PMrender::TRenderer * getRenderer();
    Plume * getPlume();
    Skirts * getSkirts();
    AirLayers * getAirLayers();
    // BrushPalette * getPalette();

    /// getter and setter for brush radii
    // float getRadius();
    // void setRadius(float rad);

    /// setter for paint or view modes
    void setMode(ControlMode mode);

    /// getter, setter, refresher for overlay texture being displayed
    void refreshOverlay();
    void setOverlay(TypeMapType purpose);
    TypeMapType getOverlay();

    /**
     * @brief bandCanopyHeightTexture   Recolour the canopy height texture according to a band of min and max tree heights
     * @param mint  Minimum tree height (below which heights are coloured black)
     * @param maxt  Maximum tree height (above which heights are coloured red)
     */
    void bandCanopyHeightTexture(float mint, float maxt);

    /**
     * Load scene attributes that are located in the directory specified
     * @param dirprefix     directory path and file name prefix combined for loading a scene
     */
    void loadScene(int numframes, int start);
    void loadScene(std::string dirprefix, int numframes, int start);

     /**
      * Save scene attributes to the directory specified
      * @param dirprefix     directory path and file name prefix combined for saving a scene, directory is assumed to exist
      */
     void saveScene(std::string dirprefix);

     /**
     * @brief readPaintMaps Input a text file encoding the paint texture layers for type and slope.
     * @param paintfile paint file name
     * @return true if the maps are read successfully, false otherwise
     */
    bool readPaintMaps(std::string paintfile);

    /**
     * @brief writePaintMaps Output text file encoding the paint texture layers for type and slope.
     * @param paintfile paint file name
     * @return true if the maps are written successfully, false otherwise
     */
    bool writePaintMaps(std::string paintfile);

    /// Add an extra scene with placeholder view, terrain and typemap onto the end of the scene list
    void addScene();

    /// change the scene being displayed
    void setScene(int s);

    /// Prepare decal texture
    void loadDecals();
    
    /// Bind the Airlayer texture to its plane
    void bindAirLayer();

    /// Load from file to appropriate TypeMap depending on purpose
    int loadTypeMap(MapFloat * map, TypeMapType purpose);

    /// Respond to key press events
    void keyPressEvent(QKeyEvent *event);

    /// set scaling value for all terrains
    void setScales(float sc);
    
    /**
      * @brief genAirLayerShape     Create geometry for planes as part of the atmosphere
      */
    void genAirLayerShape();
    
    /**
      * @brief genSkirtShape     Create geometry for spheres as part of the skirts
      */
    void genSkirtShape();
    
    /**
      * @brief genPlumeShape     Create geometry for spheres as part of the plume
      */
    void genPlumeShape();
    
    /**
      * @brief genTreeShape     Create geometry for trees with cone top
     * @param trunkheight  proportion of height devoted to bare trunk
     * @param trunkradius  radius of main trunk
      */
    void genTreeShape(float trunkheight, float trunkradius);
    
    /**
      * @brief genTreeShape     Create geometry for shrubs with sphere top
     * @param trunkheight  proportion of height devoted to bare trunk
     * @param trunkradius  radius of main trunk
      */
    void genShrubShape(float trunkheight, float trunkradius);
    
    /**
     * @brief drawPlants    Prepare plants for drawing
     * @param drawParams    Data for positioning instances
     */
    void drawPlants(std::vector<ShapeDrawData> &drawParams);
    
    /**
     * @brief drawAirLayers     Prepare volcano air layers for drawing
     * @param drawParams    Data for positioning instances
     */
    void drawAirLayers(std::vector<ShapeDrawData> &drawParams);
    
    /**
     * @brief drawPlume     Prepare volcano plume for drawing
     * @param drawParams    Data for positioning instances
     */
    void drawPlume(std::vector<ShapeDrawData> &drawParams);
    
    /**
     * @brief drawSkirts    Prepare volcano skirts for drawing
     * @param drawParams    Data for positioning instances
     */
    void drawSkirts(std::vector<ShapeDrawData> &drawParams);
    
    void setSkirtsVisible(bool show){ skirtsvisibility = show; }
    bool getSkirtsVisible(){ return skirtsvisibility; }
    
    void setPlumeVisible(bool show){ plumevisibility = show; }
    bool getPlumeVisible(){ return plumevisibility; }
    
    void setPlumeCharge(bool show){ plumecharge = show; }
    bool getPlumeCharge(){ return plumecharge; }
    
    /**
     * @brief positionPlantsToDisplay   in forest and shrub regions determine positions for geometry placement using a jittered grid
     * @param treespacing   Determines resolution of tree grid cell as treespacing x treespacing pixels
     * @param shrubspacing  Determines resolution of shrub grid cell as treespacing x treespacing pixels
     */
    void positionPlantsToDisplay(int treespacing, int shrubspacing);

signals:
    void signalRepaintAllGL();
    
public slots:
    void animUpdate(); // animation step for change of focus
    void rotateUpdate(); // animation step for rotating around terrain center

    std::string get_dirprefix();

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

    void mousePressEvent(QMouseEvent *event);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent * wheel);

private:

     QGLFormat glformat; //< format for OpenGL
    // scene control
    uts::vector<Scene *> scenes;
    std::string datadir;
    Window * mainparent;

    int currscene;
    bool dbloaded; // set to true once the user has opened a database

    // render variables
    PMrender::TRenderer * renderer;
    bool decalsbound;
    GLuint decalTexture;

    // gui variables
    bool viewing;
    bool viewlock;
    bool focuschange;
    bool focusviz;
    bool timeron;
    bool active; //< scene only rendered if this is true
    bool plumevisibility;
    bool skirtsvisibility;
    bool plumecharge;
    float scf;
    ControlMode cmode;

    QPoint lastPos;
    QColor qtWhite;
    QTimer * atimer, * rtimer; // timers to control different types of animation
    QLabel * vizpopup;  //< for debug visualisation

    // brush variables
    /*
    BrushPalette * palette;
    BrushCursor brushcursor;
    BrushPaint brush;
    */
    
    // plant geometry
    std::vector<Plant> trees;   ///< placement of trees to cover forest areas
    std::vector<Plant> shrubs;  ///< placement of shrubs to cover shrub areas
    NoiseField * nfield1, * nfield2; //< coherent noise
    Shape skirtshape;
    Shape plumeshape;
    Shape treeshape;
    Shape shrubshape;
    Shape airlayershape;

    /**
     * @brief pickInfo  write information about a terrain cell to the console
     * @param x         x-coord on terrain grid
     * @param y         y-coord on terrain grid
     */
    void pickInfo(int x, int y);
};

#endif
