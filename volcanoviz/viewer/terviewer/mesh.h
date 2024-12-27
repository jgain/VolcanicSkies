/**
 * @file
 *
 * Data structure representing a triangle mesh in 3D space.
 */

#ifndef _MESH
#define _MESH

#include <vector>
#include <stdio.h>
#include <iostream>
#include "trenderer.h"

using namespace std;

const int sphperdim = 20;

/**
 * A triangle in 3D space, with 3 indices into a vertex list and an outward facing normal. Triangle winding is counterclockwise.
 */
struct Triangle
{
    int v[3];   ///< index into the vertex list for triangle vertices
    Vector n;   ///< outward facing unit normal to the triangle
};

/**
 * An edge in 3D space, with 2 indices into a vertex list. The order of the vertices may have significance.
 */
struct Edge
{
    int v[2];   ///< indices into the vertex list for edge endpoints
};

/**
 * A triangle mesh in 3D space. Ideally this should represent a closed 2-manifold but there are validity tests to ensure this.
 */
class Mesh
{
private:
    std::vector<vpPoint> verts; ///< vertices of the tesselation structure
    std::vector<Vector> norms;  ///< per vertex normals
    std::vector<Triangle> tris; ///< triangles that join to make up the mesh
    GLfloat * col;              ///< (r,g,b,a) colour
    float scale;                ///< scaling factor
    Vector trx;                 ///< translation
    float xrot, yrot, zrot;     ///< rotation angles about x, y, and z axes

    /**
     * Search list of vertices to find matching point
     * @param pnt       point to search for in vertex list
     * @param[out] idx  index of point in list if found, otherwise -1
     * @retval true  if the point is found in the vertex list,
     * @retval false otherwise
     */
    bool findVert(vpPoint pnt, int &idx);

    /**
     * Construct a hash key based on a 3D point
     * @param pnt   point to convert to key
     * @param bbox  bounding box enclosing all mesh vertices
     * @retval hash key
     */
    long hashVert(vpPoint pnt, BoundBox bbox);

    /**
     * Construct a hash key based on the indices of an edge
     * @param v0    first endpoint index
     * @param v1    second endpoint index
     * @retval hash key
     */
    long hashEdge(int v0, int v1);

    /// Connect triangles together by merging duplicate vertices
    void mergeVerts();

    /// Generate vertex normals by averaging normals of the surrounding faces
    void deriveVertNorms();

    /// Generate face normals from triangle vertex positions
    void deriveFaceNorms();

    /**
     * Composite rotations, translation and scaling into a single transformation matrix
     * @param tfm   composited transformation matrix
     */
    void buildTransform(glm::mat4x4 &tfm);

public:

    Shape geometry;         ///< renderable version of mesh

    /// Default constructor
    Mesh();

    /// Destructor
    virtual ~Mesh();

    /// Remove all vertices and triangles, resetting the structure
    void clear();

    /// Test whether mesh is empty of any geometry (true if empty, false otherwise)
    bool empty(){ return verts.empty(); }

    /// Setter for scale
    void setScale(float scf){ scale = scf; }

    /// Getter for scale
    float getScale(){ return scale; }

    /// Setter for translation
    void setTranslation(Vector tvec){ trx = tvec; }

    /// Getter for translation
    Vector getTranslation(){ return trx; }

    /// Setter for rotation angles
    void setRotations(float ax, float ay, float az){ xrot = ax; yrot = ay; zrot = az; }

    /// Getter for rotation angles
    void getRotations(float &ax, float &ay, float &az){ ax = xrot; ay = yrot; az = zrot; }

    /// Setter for colour
    void setColour(GLfloat * setcol){ col = setcol; }

    /// Getter for number of vertices
    int getNumVerts(){ return (int) verts.size(); }

    /// Getter for number of faces
    int getNumFaces(){ return (int) tris.size(); }

    /**
     * Generate and bind triangle mesh geometry for OpenGL rendering
     * @param view      current view parameters
     * @param[out] sdd  openGL parameters required to draw this geometry
     * @retval @c true  if buffers are bound successfully, in which case sdd is valid
     * @retval @c false otherwise
     */
    // bool bindGeometry(View * view, ShapeDrawData &sdd);

    /**
     * Generate triangle mesh geometry for OpenGL rendering
     * @param[out] geom triangle-mesh geometry packed for OpenGL
     */
    void genGeometry(Shape * geom);

    /**
     * Scale geometry to fit bounding cube centered at origin
     * @param sidelen   length of one side of the bounding cube
     */
    void boxFit(float sidelen);

    /**
     * Read in triangle mesh from STL format binary file
     * @param filename  name of file to load (STL format)
     * @retval true  if load succeeds,
     * @retval false otherwise.
     */
    bool readSTL(string filename);

    /**
     * Write triangle mesh to STL format binary file
     * @param filename  name of file to save (STL format)
     * @retval true  if save succeeds,
     * @retval false otherwise.
     */
    bool writeSTL(string filename);
};

#endif
