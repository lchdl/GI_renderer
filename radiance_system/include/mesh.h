/* implements a basic model loading/saving system */

#pragma once

#include "basedefs.h"
#include "graphics.h"
#include "convex.h"
#include "basemath.h"
#include "geometry.h"
#include "stb_image/stb_image_load.h"
#include "stb_image/stb_image_write.h"
#include "random.h"

class Mesh {
protected:
    /* raw data storage */
    Array<VEC3> p; /* vertex position array */
    Array<VEC3> n; /* vertex normal array */
    Array<VEC3> t; /* vertex texture coordinate array */
                   /* for 3d tex coord, format is (u, v, w) */
                   /* for 2d tex coord, format is (u, v, 0) */
    Array<INT3> f; /* face indices (each face is a triangle) */
protected:
    /* for OpenGL */
    GLenum usage;                     /* static/stream/dynamic draw */
    GLVertexBuffer_XYZ_RGB_UV vbuf;     /* GPU vertex buffer handler */
    Array<GL_VERTEX_XYZ_RGB_UV> packedVertices; /* raw vertex buffer data stored in CPU memory */
protected:
    void packVertices(); /* assemble vertex data for uploading to GPU */
    
    enum OBJFaceFormat {
        InvalidFace,
        V,               /* vertex position only:     "f   1     2     3    " */
        V_Vt,            /* position and texcoord:    "f   3/1   4/2   5/3  " */
        V_Vn,            /* position and normal:      "f   1//2  7//8  9//3 " */
        V_Vt_Vn          /* pos, texcoord and normal: "f   1/2/3 5/6/2 3/4/5" */
    };

    enum OBJTexCoordFormat {
        InvalidTex,
        UV,              /* 2D uv coordinate  */
        UVW              /* 3D uvw coordinate */
    };

    OBJFaceFormat getOBJFaceFormat(const char* file);

public:
    Mesh();
    Mesh(const Mesh& that);
    virtual ~Mesh();

    Mesh& operator=(const Mesh& that);

    /* load Wavefront OBJ file format as mesh */
    bool loadOBJ(const char* file, GLenum usage = GL_STATIC_DRAW);
    /* assemble mesh from raw data manually */
    void assemble(
        const Array<VEC3>& p, /* positions */
        const Array<VEC3>& n, /* normals */
        const Array<VEC3>& t, /* texture coordinates */
        const Array<INT3>& f, /* face indices (v/t/n) */
        GLenum usage          /* static/stream/dynamic */
    );
    /* check if a mesh contains actual data that can be drawn onto screen */
    bool isNull() const;
    /* invoke OpenGL draw call */
    void draw();
    /* unload all resources allocated for this mesh object */
    void unload();
    /* upload data to GPU */
    void upload();
    /* retrieve member data */
    const Array<VEC3>& getPositionArray() const;
    VEC3 getVertex(const int& index) const;
    bool setVertex(const int& index, const VEC3& v);
    const Array<VEC3>& getNormalArray() const;
    const Array<VEC3>& getTexCoordArray() const;
    const Array<INT3>& getFaceIndexArray() const;
    /* caluclate CPU memory consumption */
    int bytes(); 
    /* export obj model */
    bool saveOBJ(const char* file);

    /* simple mesh operations */
    /* rotate mesh around origin (0,0,0) */
    void rotate(QUAT& q);
    /* move mesh */
    void translate(VEC3& offset);
    
    triangle getTriangle(const int& ind) const;
    triangleEx getTriangleEx(const int& ind) const;
    aabb getBoundingBox() const;

    /* mesh sampler */
    Array<VEC3> sampleMesh(const int& count) const;
    bool sampleMesh(const int& count, const char* savefile) const;
};

class MeshMaker {    
public:
    /* making simple primitives */
    static Mesh box(const VEC3& size);

    /* stitch multiple meshes into a single object */
    static Mesh stitchMeshes(const Array<Mesh*> & meshes);
};

/* make convex hull from a mesh, useful for rigid body dynamics */
Mesh buildConvexMesh(Mesh& m);

/* mesh that can be morphed using interpolations between multiple morph targets */
class MorphMesh {
protected:
    Array<Mesh> meshArray; /* array that stores multiple morph targets       */
    Mesh meshInterp;       /* interpolated mesh using multiple morph targets */
public:
    MorphMesh();
    MorphMesh(const MorphMesh& that);
    virtual ~MorphMesh();

    MorphMesh& operator=(const MorphMesh& that);
    
    /* 

    Desription
    ----------    
    Load a mesh file as morph target, returns an integer   
    representing the handle of loaded morph target (>0).    
    "0" will be returned if any error occurs when loading. 

    Note 
    ----------
    Each morph target should share the same topology, 
    otherwise error will occur when interpolating. 
    The first loaded morph target is considered as the base 
    mesh. If no interpolation has done, the first loaded   
    mesh is used.                                          

    */
    unsigned int loadMorphTarget(const char* file);

    /* 

    Desription
    ----------
    Calculate interpolated (morphed) mesh using multiple morph targets.

    Parameters
    ----------
    n: number of morph targets used in interpolation  
    meshHandles: mesh handles returned from loadMorphTarget(...) 
                 indicating which meshes are used for morphing. 
    mixWeights: mixing weights of each mesh. Value    
                range is not limited to [0,1] for flexbility.
    
    */
    bool interpMeshes(int n, unsigned int* meshHandles, REAL* mixWeights);

    /* unload all morph targets and interpolated mesh */
    void unload();

    /* draw interpolated mesh */
    void draw();
};

/* static mesh with BST hierarchical structure */
#define BSTMESH_DBG_LOGBUFFER_SIZE   1024
class BSTMesh : public Mesh, public Noncopyable{
protected:
    enum BSTMeshNodeType {
        Invalid,
        NonLeaf,
        OrdinaryLeaf, /* depth <= max_depth && tris_count <= max_count */
        AtomLeaf,     /* leaves that cannot be split effectively */
        DeepLeaf,     /* depth == max_depth && tris_count > max_count */
    };
    struct BSTMeshNode : public Noncopyable {
    public:
        aabb box;                   /* axis aligned bounding box */
        BSTMeshNode* left, *right;  /* left and right node */
        int depth;                  /* node depth */
        int n;                      /* number of triangles in this node (leaf) */
        int* idxs;                  /* the first vertex indices of all the triangles in this node */
        BSTMeshNodeType type;       /* node type (used in tree structure analysis) */
    public:
        BSTMeshNode();
        virtual ~BSTMeshNode();
        bool isLeaf() const;
    };
    struct BSTMeshStatistics {
        int nTotalNodes;
        int nInvalidNodes;
        int nNonLeaf;
        int nOrdinaryLeaf;
        int nAtomLeaf;
        int nDeepLeaf;
        Array<int> leafTris; /* number of triangles in all leaf nodes */
        int nLeafTrisSum;    /* sum of all triangles in all leaf nodes */
        int nMaxLeafTris;    /* = max(leafTris) */
        REAL overlapRatio;   /* nLeafTrisSum / num_of_tris, smaller value indicates better performance */
        int nBytes;

        BSTMeshStatistics();
        void clear();
    };
    struct BSTMeshRenderBuffer {
        int w, h;
        VEC3* pixels;
    };
public:
    struct BSTMeshBuildSettings {
        int maxDepth;       /* max node depth */
        int splitThreshold; /* split node if triangle count is larger than this value */
        int sahResolution;  /* higher resolution will result in more refined tree structure, */
                            /* but will increase build time dramatically. */
        BSTMeshBuildSettings();
        bool isValid();
    };
    static BSTMeshBuildSettings defaultBuildSettings();
protected:
    BSTMeshBuildSettings bstSettings;
    BSTMeshNode*         bstRoot;     /* BST root node */
    BSTMeshStatistics    bstStat;
    GLShader               bstShader;   /* shader for drawing the wireframe */
    GLVertexBuffer_XYZ     bstVB;

protected:
    void buildBST_splitBox(const aabb* box, const char axis, const REAL& value, aabb* left, aabb* right);
    void buildBST_sahCost(const Array<int>& idxs,
        const aabb& left, const aabb& right,
        const char& axis, const REAL& value,
        REAL * cost);
    void buildBST_recur(BSTMeshNode* node, const Array<int>& idxs);
    void buildBST_boxFit(const Array<int>& idxs, VEC3* bmin, VEC3* bmax);
    void buildBST(const BSTMeshBuildSettings& settings);
    void clearBST_recur(BSTMeshNode* node);
    void clearBST();

public:
    BSTMesh();
    virtual ~BSTMesh();

    /* load Wavefront OBJ file format as mesh */
    bool loadOBJ(const char* file, const BSTMeshBuildSettings& settings, GLenum usage = GL_STATIC_DRAW);
    /* build mesh manually */
    void assemble(
        const Array<VEC3>& p, /* positions */
        const Array<VEC3>& n, /* normals */
        const Array<VEC3>& t, /* texture coordinates */
        const Array<INT3>& f, /* face indices (v/t/n) */
        const BSTMeshBuildSettings& settings, /* BST build settings */
        GLenum usage          /* static/stream/dynamic */
    );
    /* unload all resources allocated for this mesh object */
    void unload();
    /* simple mesh operations */
    /* rotate mesh around origin (0,0,0) */
    void rotate(QUAT& q);
    /* move mesh */
    void translate(VEC3& offset);

    /* I/O */
protected:
    void saveNode_recur(BSTMeshNode* node, FILE* fp);
    void loadNode_recur(BSTMeshNode** node, FILE* fp);
public:
    bool saveAsset(const char* file);
    bool loadAsset(const char* file);

    /* advanced */
protected:
    void bstAdvInit();
    void bstAdvFree();

    void drawNodeBoundary_recur(GLCamera* camera, BSTMeshNode * node, int minDepth, int maxDepth);
    void doStatistics_recur(BSTMeshNode * node);
    BYTE doRT_byteClamp(const REAL& value); /* clamp REAL range to 0~255 */

public:
    void drawNodeBoundary(GLCamera* camera, const VEC3& color, int minDepth, int maxDepth);
    void doStatistics(char** log); /* please pass a NULL pointer to retrieve log!!! */
    void doRaytracing(const VEC3& pos, const VEC3& look, const VEC3& up, const REAL& fov,
        const int& w, const int& h, const char* image);
    void doRaytracing(const ray& r, ray_hit& rtResult, const REAL& searchDist);
    void doRaytracing(const ray& r, Array<ray_hit> & rtResults, const REAL& searchDist);
    void doSphereCrossing(const sphere& s, Array<int>& idxs);
};
