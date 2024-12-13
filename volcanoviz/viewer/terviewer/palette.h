/**
 * @file
 */

// palette.h: ecosystem painting controls to apply plant distributions onto a terrain
// author: James Gain
// date: 4 July 2016

#ifndef _constraint_h
#define _constraint_h

#include "stroke.h"
#include "typemap.h"
#include "shape.h"
#include <common/debug_vector.h>
#include <common/debug_list.h>
#include <QWidget>
#include <QPushButton>
#include <QImage>
//#include <common/map.h>

const float manipradius = 75.0f;
const float manipheight = 750.0f;
const float armradius = manipradius / 2.5f;
const float tolzero = 0.01f;

class ConditionsMap;
class GLWidget;

enum class ControlMode
{
    VIEW,   // free viewing of terrain
    PLACEHERD,  // paint the herd starting locatoin
    PAINTBOIDS, // painting boid behaviour
    PAINTTYPE,  // painting terrain type
    PAINTSLOPE, // painting terrain slope categories
    CMEND
};

const std::array<BrushType, 21> all_brushlandtypes = {BrushType::FREE, BrushType::ROCK, BrushType::TRAIL, BrushType::GRASS, BrushType::SHRUBS, BrushType::FOREST, BrushType::WATER, BrushType::SNOW, BrushType::ICE, BrushType::OBSTACLE, BrushType::GENTLE, BrushType::STEEP, BrushType::CLIFF,
    BrushType::BOIDS1, BrushType::BOIDS2, BrushType::BOIDS3, BrushType::BOIDS4, BrushType::BOIDS5,
    BrushType::BOIDS6, BrushType::HERDADD, BrushType::HERDDEL
}; // to allow iteration over the brushtypes

#define PALETTE_ENTRIES 21

class BrushPalette : public QWidget
{
    Q_OBJECT

public:

    BrushPalette(TypeMap * typemap, QWidget *parent = 0);

    ~BrushPalette(){}

    //QSize minimumSizeHint() const;
    QSize sizeHint() const;

    /// Obtain the current active brush type from the palette
    BrushType getDrawType(){ return typeSel[currSel]; }

    /// move the currently active icon to the correct palette entry
    void setActivePalette();

    /// Set the currently active brush type in the palette
    void setDrawType(BrushType btype);

    /// Enable or Disable palette brush icon
    void setDisabled(BrushType btype);
    void setEnabled(BrushType btype);
    
    /// Report whether a palette entry is present or not depending on the particular species
    bool inSpeciesPalette(int entry);
    
    /// Return the initially active palette
    BrushType firstSpeciesPalette();

public slots:

    /// palette entry button press
    void typeSelect();

    /**
     * @brief directSelect Directly set a particular brush to active
     * @param brushtype The brush to be activated
     */
    void directSelect(BrushType brushtype);

    /**
     * @brief switchMode Display correct palette entries depending on painting mode
     * @param cmode New painting mode to activate
     */
    void switchMode(ControlMode cmode);

private:
    GLWidget * glparent;
    TypeMap * tmap;
    QImage unassignedImg, activeImg, forestImg, shrubImg, waterImg, grassImg, iceImg, snowImg, obstacleImg;
    QPushButton * selector[PALETTE_ENTRIES];
    QPushButton * resetor[PALETTE_ENTRIES];
    BrushType typeSel[PALETTE_ENTRIES];
    int currSel;
};


/// Widget that displays a double torus over the landscape to indicate bounds of a paint operation
class BrushCursor
{
private:
    vpPoint pos;        ///< cursor position while over terrain
    // vpPoint trackpnt;   ///< debugging
    bool active;        ///< determines whether the radius indicator should be displayed
    float radius;  ///< radial effect of painting on terrain
    float hghtoffset;  ///< offset to raise ring so that it clear of the terrain

public:

    Shape shape; ///< geometry for rendering
    // Shape debugshape;

    BrushCursor(){ active = false; hghtoffset = 0.0f; }

    /// setter for active indicator
    void setActive(bool on){ active = on; }

    /// getter for active indicator
    bool getActive(){ return active; }

    /// getter and setter for brush radii
    void setRadius(float rad);
    float getRadius(){ return radius; }

    /// setter for height offset
    void setHeightOffset(float off){ hghtoffset = off; }

    /// setter for brush colour
    void setBrushColour(GLfloat * col)
    {
        shape.setColour(col);
    }

    /**
     * Create manipulator geometry on update
     * @param view      current view state
     * @param terrain   terrain being synthesized
     * @param brushradius radius of ring
     * @param dashed    whether or not to draw in a dashed style
     */
    void genBrushRing(View * view, Terrain * terrain, float brushradius, bool dashed);

    /**
     * Update the cursor position for rendering the brush radius
     * @param view      current view state
     * @param terrain   terrain being synthesized
     * @param x, y      on-screen mouse position
     * @retval @c true if the cursor is over the terrain
     */
    bool cursorUpdate(View * view, Terrain * terrain, int x, int y, int port[4]);

    /**
     * Get terrain type corresponding to current mouse screen coordinates by picking
     * @param view      current view state
     * @param terrain   terrain being synthesized
     * @param tmap      terrain type map
     * @param x, y      on-screen mouse position
     */
    // int pickType(View * view, Terrain * terrain, TypeMap * tmap, int x, int y);
};

class BrushPaint
{
private:
    Terrain * terrain;      ///< to access the terrain
    Region coverage;        ///< bounding box for stroke coverage
    float radius;           ///< radial effect of stroke on terrain
    BrushType brushtype;    ///< action of brush on the map
    vpPoint currpnt;        ///< most recent mouse position on the terrain
    vpPoint prevpnt;        ///< previous painted mouse position on the terrain
    bool drawing;           ///< is drawing active
    BoundRect bnd;          ///< bounding box for update to paint map

    /**
     * Write the brush stroke col to the Ecosys Paint Map
     * This is also where ecosys updates will happen in due course
     * @param pmap          type map holding all paint colours
     * @param radius        radial effect of the brush
     */
    void paintMap(TypeMap * pmap, float radius);

public:

    BrushPaint(){ brushtype = BrushType::GRASS; drawing = false; }

    /**
     * Create manipulator geometry on update
     * @param view  current view state
     */
    void genManipulator(View * view);

    /**
     * Constructor
     * @param ter        terrain being synthesized
     * @param btype      type associated with brush
     */
    BrushPaint(Terrain * ter, BrushType btype);

    ~BrushPaint(){}

    /// setter for brush type
    void setBrushType(BrushType type){ brushtype = type; }

    /**
     * Add a point in screen coordinates to the current stroke and update the paint map accordingly
     * @param view      current view state
     * @param pmap      paint map
     * @param x, y      on-screen mouse position
     * @param radius    radial effect of the brush
     * @param port      viewport parameters
     */
    void addMousePnt(View * view, TypeMap * pmap, int x, int y, float radius, int port[4]);

    /**
     * @brief startStroke Start the current brush stroke. Called on mouse down.
     */
    void startStroke();

    /**
     * Complete the current brush stroke. Called on mouse up.
     */
    void finStroke();
};

#endif
