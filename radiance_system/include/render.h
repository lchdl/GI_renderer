/* define scenes, images, and all kinds of renderers */
#pragma once

#include "basedefs.h"
#include "basemath.h"
#include "stb_image/stb_image_load.h"
#include "stb_image/stb_image_write.h"
#include "mesh.h"
#include "geometry.h"
#include "random.h"
#include "graphics.h"
#include "timer.h"

typedef void(*pfnColorFunction)(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a);
void vflipBuffer(BYTE * p, int stride, int h);
BYTE byteClamp(const REAL& value);
void albedoToRGBA(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a);
void normalToRGBA(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a);
VEC3 computeReflectionDir(VEC3 n, VEC3 d);
VEC3 computeRefractionDir(VEC3 d, VEC3 n, REAL ior_in, REAL ior_out);
REAL clampComp(const REAL& v, const REAL& compMin = ZERO, const REAL& compMax = ONE);
VEC3 clampComp(const VEC3& v, const REAL& compMin = ZERO, const REAL& compMax = ONE);
VEC4 clampComp(const VEC4& v, const REAL& compMin = ZERO, const REAL& compMax = ONE);
VEC3 ARGBToRGB(const VEC4& v);
int ipow(int base, int power); /* utility function to compute base^power */
REAL computeSmoothness(const VEC3& n1, const VEC3& n2);

/* normal mapping */
void computeNormalMappingTBN(
    const VEC3& p1, const VEC3& p2, const VEC3& p3,
    const VEC2& t1, const VEC2& t2, const VEC2& t3,
    VEC3& T, VEC3& B, VEC3& N);
void computeNormalMappingTBN(const triangleEx& tri, VEC3& T, VEC3& B, VEC3& N);
void normalMappingTransformToWorld(const VEC3& T, const VEC3& B, const VEC3& N,
    const VEC3& localNormal, VEC3& worldNormal);

/* utility function */
VEC3 rgbToVec3(BYTE r, BYTE g, BYTE b);
VEC4 rgbToVec4(BYTE r, BYTE g, BYTE b);
VEC4 argbToVec4(BYTE a, BYTE r, BYTE g, BYTE b);

/* here i list all base classes */
struct Material;
class  MapSampler;
class  Entity;
class  Scene;
class  InteractiveRenderer;

/* all resources used in rendering */
struct SceneResources : public Noncopyable {
    Array<MapSampler*> maps;
    Array<Material*> materials;
    Array<Entity*> entities;
    Array<int> emissiveEntityIds;
};

/* used in raytrace */
class ray_query {
protected:
    Array<ray_hit> hits;
public:
    void addHit(ray_hit& hit);
    void sort(); /* sort records based on distance */
    int numHits() const;
    ray_hit& operator [](const int & idx);
    void clear();
    void removeFirstN(const int& N); /* remove first N hits */
};

/* define map sampler texturing */
class MapSampler {
public:
    virtual VEC4 sample(VEC3 uvw) = 0;
    virtual const char* getClassName() const = 0;
};

/* material definition */
/* all material classes should not contain any pointer, */
/* they should just be a pure data structure */
enum MaterialType {
    Material_Emissive = 0x01,
    Material_Diffuse = 0x02,
    Material_Reflect = 0x04,
    Material_Refract = 0x08,
    Material_Volumetric = 0x10,
};
struct Material {

protected:
    /*
    EMISSIVE: (mutually exclusive with all other types)
    NOTE: if emissiveAmount > 0, then DIFFUSE, REFLECT, REFRACT, and VOLUMETRIC
          settings are all ignored.
    */
    VEC3 emissiveIntensity; /* RGB (0 ~ +inf) */
    int  emissiveMapId;
    int  emissiveSamplesNum;
    REAL emissiveFalloff; /* 0 ~ +inf */

    /*
    DIFFUSE / REFLECT:
    */
    VEC3 diffuseColor; /* RGB (0~1) */
    int  diffuseMapId; /* resources.maps[diffuseId] >>> MapSampler* */
    int  normalMapId;
    REAL diffuseRoughness;
    VEC3 reflectColor;
    int  reflectMapId; /* resources.maps[diffuseId] >>> MapSampler* */
    REAL reflectRoughness;
    bool fresnelReflect;

    bool thinObject;
    VEC3 thinObjectFalloff;

    /*
    REFRACT:
    */
    VEC3 refractColor;
    REAL refractIndex; /* refractive index, IOR */
    REAL refractRoughness;

    /*
    VOLUMETRIC / SSS:
    */
    VEC3 volumeColor;
    REAL volumeFalloff;

public:

    Material();
    virtual ~Material();

    unsigned int getMaterialType() const;

    VEC4 getDiffuseColor(VEC3 uvw, const Scene* scene) const;
    VEC3 getLocalNormal(VEC3 uvw, const Scene* scene) const;
    int  getNormalMapId() const;
    VEC4 getEmissiveIntensity(VEC3 uvw, const Scene* scene) const;
    REAL getEmissiveFalloff() const;
    int  getNumEmissiveSamples() const;
    REAL getRefractionIndex() const;

    void setEmissive(VEC3 emissiveColor, int emissiveMapId, REAL emissiveFalloff, int emissiveSamples = 1000);
    void setDiffuse(VEC3 diffuseColor, int diffuseMapId, REAL diffuseRoughness, int normalMapId, 
        bool thinObject = false, VEC3 thinObjectFalloff = VEC3::ones());
    void setReflect(VEC3 reflectColor, int reflectMapId, REAL reflectRoughness, bool fresnelReflect);
    void setRefract(VEC3 refractColor, REAL refractIndex, REAL refractRoughness);
    void setVolumetric(VEC3 volumeColor, REAL volumeFalloff);

    bool hasProperty(unsigned int materialTypeMask);
    bool isThinObject() const;
    VEC3 getThinFalloff() const;
};

class Transform {
protected:
    /* transform order: scale -> rotate -> translate */
    REAL scale;     /* proportional scaling */
    QUAT rotation;
    VEC3 translation;
public:
    Transform();
    Transform(const REAL& scale, const QUAT& rotation, const VEC3& translation);
    virtual ~Transform();
    Transform operator =(const Transform&);

    void setTransform(const REAL& scale, const QUAT& rotation, const VEC3& translation);
    void setTransform(const Transform& xform);
    VEC3        transformPointToLocal(const VEC3& p) const;
    Array<VEC3> transformPointToLocal(const Array<VEC3>& p) const;
    VEC3        transformVectorToLocal(const VEC3& v) const;
    Array<VEC3> transformVectorToLocal(const Array<VEC3>& v) const;
    VEC3        transformPointToWorld(const VEC3& p) const;
    Array<VEC3> transformPointToWorld(const Array<VEC3>& p) const;
    VEC3        transformVectorToWorld(const VEC3& v) const;
    Array<VEC3> transformVectorToWorld(const Array<VEC3>& v) const;
    aabb transformBoundingBoxToWorld(const aabb& bbox) const;

    REAL getScale() const;
    QUAT getRotation() const;
    VEC3 getTranslation() const;

};

/* define interface for all renderable objects */
class Entity : public Noncopyable {
    /*
    Since we have proxy mechanics, we do not need to deep copy an entity.
    For safety reason I disabled the copy constructor and assignment operator.
    */
protected:
    Array<VEC3> lightSamples;
    Transform xform;             /* local to world transform */
    aabb bboxWorld;              /* world bounding box */
    int materialId;
    bool visible;

protected:
    void updateWorldBoundingBox();

public:
    Entity();
    virtual ~Entity();

    virtual const char* getClassName() const = 0;

    virtual Array<ray_hit> rayQuery(const ray* r, const REAL& searchDist) = 0;
    virtual Array<ray_hit> rayQuery(const ray* r, const Transform* proxyXform, const REAL& searchDist) = 0;

    void useMaterial(const int& materialId);
    int  getMaterialId() const;

    /* for normal mapping */
    virtual VEC3 normalQuery(const int& triangleId, const VEC3& uvw, const Scene* scene) const = 0;
    virtual VEC3 normalQuery(const int& proxyMaterialId, const Transform* proxyXform, 
        const int& triangleId, const VEC3& uvw, const Scene* scene) const = 0;

    /* define the bounding box of each entity, used in BVH construction */
    aabb getWorldBoundingBox();
    virtual aabb getLocalBoundingBox() const = 0;

    /* visibility */
    bool isVisible() const;
    void setVisible(bool state = true);

    /* emissive property */
    virtual void buildLightSamples(const int& count) = 0;
    VEC3         getLocalLightSample(XorwowRNG& rng) const; /* get a single light sample position in local space */
    virtual VEC3 getWorldLightSample(XorwowRNG& rng) const = 0;
    virtual VEC3 getWorldLightSample(XorwowRNG& rng, const Transform* proxyXform) const = 0;

    /* transform */
    void setTransform(const Transform& xform);
    void setTransform(const REAL& scale, const QUAT& rotation, const VEC3& translation);
    Transform getTransform() const;
};
class EntityProxy : public Entity {
protected:
    Entity* proxyTarget;

public:
    EntityProxy();
    virtual ~EntityProxy();

    virtual const char* getClassName() const;

    void setProxyTarget(Entity* entity);

    virtual VEC3 normalQuery(const int& triangleId, const VEC3& uvw, const Scene* scene) const;
    virtual VEC3 normalQuery(const int& proxyMaterialId, const Transform* proxyXform, 
        const int& triangleId, const VEC3& uvw, const Scene* scene) const;

    virtual void buildLightSamples(const int& count);
    VEC3         getLocalLightSample(XorwowRNG& rng) const;
    virtual VEC3 getWorldLightSample(XorwowRNG& rng) const;
    virtual VEC3 getWorldLightSample(XorwowRNG& rng, const Transform* proxyXform) const;

    virtual aabb getLocalBoundingBox() const;
    virtual Array<ray_hit> rayQuery(const ray* r, const REAL& searchDist);
    virtual Array<ray_hit> rayQuery(const ray* r, const Transform* proxyXform, const REAL& searchDist);

};
class EntityStack : protected Stack<int>
{
public:

    EntityStack();
    virtual ~EntityStack();

    void pushEntityId(const int& entityId);
    int popEntityId();

    REAL getCurrentIOR(const Scene* scene) const;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

class Bitmap : public MapSampler
{
public:

    Bitmap();
    virtual ~Bitmap();
    Bitmap(const Bitmap&);
    Bitmap operator=(const Bitmap&);
    bool load(const char* file, bool vflip = true);
    bool save(const char* file, bool vflip = false);
    void create(int w, int h);
    void unload();
    int getWidth() const;
    int getHeight() const;

protected:
    int toOffset(int x, int y);
    VEC4 toVec4(unsigned int x);

public:
    /* image sampling */
    virtual VEC4 sample(VEC3 uvw);
    virtual const char* getClassName() const;

    static Bitmap stitchBitmaps(const Array<String>& files);

protected:
    int         w, h, c; /* width, height, channels */
    void*       data;

};

class MeshEntity : public Entity, public BSTMesh {
public:
    MeshEntity();
    virtual ~MeshEntity();

    virtual const char* getClassName() const;

    bool loadAsset(const char* file);

    virtual void buildLightSamples(const int& count);
    virtual VEC3 getWorldLightSample(XorwowRNG& rng) const;
    virtual VEC3 getWorldLightSample(XorwowRNG& rng, const Transform* proxyXform) const;

    virtual VEC3 normalQuery(const int& triangleId, const VEC3& uvw, const Scene* scene) const;
    virtual VEC3 normalQuery(const int& proxyMaterialId, const Transform* proxyXform, const int& triangleId, const VEC3& uvw, const Scene* scene) const;

    virtual aabb getLocalBoundingBox() const;
    virtual Array<ray_hit> rayQuery(const ray* r, const REAL& searchDist);
    virtual Array<ray_hit> rayQuery(const ray* r, const Transform* proxyXform, const REAL& searchDist);
};

struct Camera
{
    VEC3 pos;
    VEC3 look;
    VEC3 up;
    REAL fov; /* in radians */
    REAL focal; /* focal distance */

    ray raycast(const int& w, const int& h, const int& x, const int& y, 
        bool randomSampling = false, const int& subpixels = 1, XorwowRNG* rng = NULL);
    REAL convertPixelSize(int px, int w, REAL viewDistance);
    
    void worldPositionToScreen(VEC3 p, int w, int h, int& x, int& y);
};

struct light_cache
{
    VEC3 p; /* sample point */
    VEC3 n; /* surface normal */
    VEC3 i; /* sample intensity */
};

class LightCacheSystem
{
protected:
    Array<light_cache> lightCache;

protected:

    struct LightCacheKdTreeNode
    {
        aabb box;                   /* axis aligned bounding box */
        LightCacheKdTreeNode* left, * right;  /* left and right node */
        int depth;                  /* node depth */
        int n;                      /* number of light cache in this node (leaf) */
        int* idxs;                  /* light cache index */

        LightCacheKdTreeNode();
    };

    LightCacheKdTreeNode* kdroot;
    int kdMaxDepth;
    int kdSplitThreshold;
    int kMax;
    int _buildCounter;

    void sortLightCache(int start, int count, int axis);
    
    void buildKdTree_splitBox(const aabb* box, const char axis, const REAL& value, aabb* left, aabb* right);
    void buildKdTree_boxFit(const int& start, const int& count, VEC3* bmin, VEC3* bmax);
    void buildKdTree_recur(LightCacheKdTreeNode* node, const int& start, const int& count);
    void destroyKdTree_recur(LightCacheKdTreeNode* node);

    void knn_findInitialSolution(const VEC3& query, const int& k, 
        PriorityQueue<int, REAL>& solution, LightCacheKdTreeNode** searchedLeaf) const;
    void knn_refineSolution(const sphere& sph, const int& k, const LightCacheKdTreeNode* searchedLeaf,
        PriorityQueue<int, REAL>& solution) const;

public:
    LightCacheSystem();
    virtual ~LightCacheSystem();

    void addLightCache(const light_cache& cache);
    void addLightCache(const Array<light_cache>& cache);
    void buildKdTree(const int& k, const int& maxDepth); /* build kd-tree */
    void destroyKdTree();
    bool save(const char* file) const;
    bool load(const char* file);
    Array<light_cache> knnSearch(const VEC3& query, const int& k, bool approximate = false) const;

    int numCaches() const;
    int computeCacheSize() const;
    int maxNodeCacheCount() const;
    int minNodeCacheCount() const;
    const light_cache& getCache(const int& index) const;

};

class Scene : public Noncopyable
{
protected:
    /* BVH containing the full entities (renderable objects) */
    struct BVHNode {
        aabb box;
        int depth;
        int n;          /* number of entities in this node (only used in leaf node) */
        int* entityIds; /* all entity IDs */
        BVHNode *left, *right;
    public:
        BVHNode();
        virtual ~BVHNode();
        bool isLeaf() const;
    };
public:
    struct SceneBuildSettings {
        int maxDepth;       /* max node depth */
        int splitThreshold; /* split node if entity count is larger than this value          */
        int sahResolution;  /* higher reoslution will result in more refined tree structure, */
                            /* but will increase build time dramatically. */
        SceneBuildSettings();
        bool isValid();
    };
    static SceneBuildSettings defaultBuildSettings();
protected:
    void buildBVH_splitBox(const aabb* box, const char axis, const REAL& value, 
        aabb* left, aabb* right);
    void buildBVH_sahCost(const Array<int>& idxs,
        const aabb& left, const aabb& right,
        const char& axis, const REAL& value,
        REAL * cost);
    void buildBVH_recur(BVHNode* node, const Array<int>& idxs);
    void buildBVH_boxFit(const Array<int>& idxs, VEC3* bmin, VEC3* bmax);
    void buildBVH(const SceneBuildSettings& settings);
    void clearBVH_recur(BVHNode* node);
    void clearBVH();

protected:
    SceneResources     resources;
    SceneBuildSettings bvhSettings;
    BVHNode*           bvhRoot;
    Array<Camera>      cameraPresets;
public:
    int  addMeshEntity(const char* file); /* load a static mesh asset into scene */
    int  loadBitmap(const char* file);
    int  addMaterial(const Material* material);
    int  addCameraPreset(const Camera* camera);
    bool setEntityMaterial(const int& entityId, const int& materialId);
    void setEntityTransform(const int& entityId, const Transform& xform);
    void setEntityTransform(const int& entityId, const VEC3& translation = VEC3::zeros(), const QUAT& rotation = QUAT::identity(), const REAL & scale = ONE);
    int  copyEntity(const int& entityId);
    bool buildScene(const SceneBuildSettings& settings);
    void unload();    /* unload full scene */

    int        getEntityCount() const;
    const Camera& getCameraPreset(const int& cameraId) const;
    Entity*    getEntity(const int& entityId) const;
    Material*  getEntityMaterial(const int& entityId) const;
    int        getMaterialCount() const;
    Material*  getMaterial(const int& materialId) const;
    Array<int> getEmissiveEntityIds() const;
    const SceneResources& getResources() const;

    ray_query  rayQuery(const ray* r, const REAL searchDist = REAL_MAX) const;
    ray_query  rayQueryNormalMapping(const ray* r, const REAL searchDist = REAL_MAX) const;
    ray_query  rayQueryIgnoreMaterial(const ray* r, const unsigned int& ignoredMaterial, const REAL searchDist = REAL_MAX) const;
    ray_query  rayQueryIgnoreEntity(const ray* r, const unsigned int& ignoredEntityId, const REAL searchDist = REAL_MAX) const;
    

    /* GI */
    VEC3 computeDirectLightIntensity(const VEC3& samplePoint, const VEC3& normal, const int& samplesPerLight, XorwowRNG& rng);
    VEC3 finalGather(const VEC3& samplePoint, const VEC3& normal, 
        const int& hemiSamples, const int& k, const LightCacheSystem& lightCacheSys, 
        const bool& checkVisibility,
        XorwowRNG& rng);
    bool checkVisibilityBetween(const VEC3& p1, const VEC3 p2);
    bool checkVisibilityBetween(const VEC3& p1, const VEC3 p2, const VEC3& n1, const VEC3& n2);

    Scene();
    virtual ~Scene();
};

class RenderBuffer
{
protected:
    int w, h;
    VEC4* buffer;
    BYTE* raw; /* raw RGBA data */

protected:

public:
    RenderBuffer();
    virtual ~RenderBuffer();
    RenderBuffer(const RenderBuffer&);

    void destroy();
    void create(int w, int h);
    void setPixel(int x, int y, const VEC4& value);
    void setBlock(int x, int y, int w, int h, const VEC4& value);
    VEC4 getPixel(int x, int y) const;
    int getWidth() const;
    int getHeight() const;

    RenderBuffer operator + (const RenderBuffer&);
    void operator += (const RenderBuffer&);
    RenderBuffer operator * (const REAL&);
    RenderBuffer operator / (const REAL&);
    void operator /= (const REAL&);

    RenderBuffer operator = (const RenderBuffer&);

    void saveBitmap(pfnColorFunction func, const char* file);
    void saveBitmap(pfnColorFunction func, REAL scale, const char* file);
    const BYTE* asBitmap(pfnColorFunction func);
    const BYTE* asBitmap(pfnColorFunction func, REAL scale);

    void collectRegion(int x, int y, int w, int h, Array<VEC4>& pixels, 
        bool keepShape = false) const;

    VEC4 max() const;
};
struct RenderSettings
{
    /* resolution & sampling */
    struct _Render {
        VEC4 backgroundColor;
        int  w, h;
        bool antialiasing;
    } render;
    /* GI */
    struct _GI 
    {
        int  traceDepth;         /* trace depth of GI rays */
        /* resolution of GI rays, 0 means each pixel will shoot one GI ray */
        /* into the scene, -1 means every 2x2 pixel will shoot one GI ray */
        /* -2: 4x4, -3: 8x8, ... */
        int  resolution_directLighting_min;
        int  resolution_directLighting_max;
        int  resolution_indirectLighting_min;
        int  resolution_indirectLighting_max;
        bool randomSampling;     /* random (jitter) sampling from screen */
        int  numSamplesPerLight; /* number of sample rays shoot from object surface to each light */
        /* determine how many samples will be found in a k-NN search */
        int  knnSamples_directLighting;
        int  knnSamples_indirectLighting;
        int  maxKdTreeDepth;      /* kd-tree max depth */
        int  numFinalGatherSamples;
        bool checkVisibility;     /* enable visibility check during final gathering */
    } GI;

    void validateParameters();

    static RenderSettings defaultRenderSettings();

};
enum RenderStatus {
    Render_NoTask,
    Render_Running,
    Render_Finished,
};
class RenderTask
{
    friend class InteractiveRenderer;
protected:
    GLTexture2D displayTexture;
    Timer timer;
    /* utility functions */
    unsigned int toRGBA(BYTE r, BYTE g, BYTE b);
    unsigned int toRGBA(BYTE r, BYTE g, BYTE b, BYTE a);
    void drawRect(BYTE* p, int x, int y, int w, int h, unsigned int color, int W, int H);
    void drawPoint(BYTE* p, int x, int y, unsigned int color, int W, int H);
    String formatTime(double sec);
    String formatProgress(const int& total, const int& finished);
    String formatSize(int bytes);
public:
    virtual void         interactiveSubmit(InteractiveRenderer* renderer) = 0;
    virtual RenderStatus interactiveRender(InteractiveRenderer* renderer) = 0;
    virtual void         onRenderStart(InteractiveRenderer* renderer);
    virtual void         onRenderFinish(InteractiveRenderer* renderer);
    virtual bool         checkRenderTask(String& log);
};
class InteractiveRenderer {
    /* * * * * * * * * * * * * */
    /* frontend display system */
    /* * * * * * * * * * * * * */
protected:
    GLTrueTypeFont sysFont; /* used for display system messages */
    GLTextureBlitter2D_RGBA texBlit;
    GLViewport mainViewport;
    int windowWidth, windowHeight; /* main window width/height */
    bool fatalError;
protected:
    struct LogItem {
        String msg;
        int level; /* 0: info, 1: warning, 2: error */
    };
    FixedQueue<LogItem> logs;
    void logMessage(const char* msg, const int& level = 0);
public:

    void initFrontendDisplay(int w, int h);
    void onWindowResize(int w, int h);

    void logInfo(const char* msg, bool suppressDuplicates = true);
    void logWarning(const char* msg, bool suppressDuplicates = true);
    void logError(const char* msg, bool suppressDuplicates = true);
    void updateLog(const char* msg, const int& level = 0); /* update the latest line printed on screen */
    int getWindowWidth() const;
    int getWindowHeight() const;

    /* * * * * * * * * * */
    /* rendering system  */
    /* * * * * * * * * * */
public:
    InteractiveRenderer();
    virtual ~InteractiveRenderer();
protected:
    FixedQueue<RenderTask*> renderTasks;
    RenderStatus prevStatus;
    int renderCalls;
    void display(); /* display all contents onto screen */
public:
    void addRenderTask(RenderTask* task);
    bool checkRenderTasks();
    RenderTask* getCurrentRenderTask() const;
    RenderStatus interactiveRender();
    void setErrorFlagAndQuit(const String& errMsg);

};
class RenderTask_Albedo : public RenderTask
{
protected:
    Scene* scene;
    Camera camera;
    RenderSettings settings;
    RenderBuffer renderBuffer;
    String savePath;
    struct InternalVariables {
        int scan_y;
    } states;
public:
    virtual void saveImage(const char* file);
    virtual void createRenderTask(Scene* scene, Camera* camera, RenderSettings* settings, const char* save = NULL);
    virtual RenderStatus interactiveRender(InteractiveRenderer* renderer);
    virtual void interactiveSubmit(InteractiveRenderer* renderer);
    virtual void onRenderStart(InteractiveRenderer* renderer);
    virtual void onRenderFinish(InteractiveRenderer* renderer);
    RenderBuffer getRenderBuffer() const;
};
class RenderTask_Normal : public RenderTask_Albedo
{
public:
    virtual void saveImage(const char* file);
    virtual RenderStatus interactiveRender(InteractiveRenderer* renderer);
    virtual void interactiveSubmit(InteractiveRenderer* renderer);
    virtual void onRenderStart(InteractiveRenderer* renderer);
};
class RenderTask_Depth : public RenderTask_Albedo
{
public:
    virtual void saveImage(const char* file);
    virtual RenderStatus interactiveRender(InteractiveRenderer* renderer);
    virtual void interactiveSubmit(InteractiveRenderer* renderer);
    virtual void onRenderStart(InteractiveRenderer* renderer);
};
class RenderTask_GI : public RenderTask 
{
protected:
    RenderTask_Albedo* albedoTask;
    RenderTask_Normal* normalTask;
    RenderTask_Depth*  depthTask;
    RenderBuffer       renderBuffer;
    RenderBuffer       detailBuffer;
    Scene*           scene;
    Camera           camera;
    RenderSettings   settings;
    XorwowRNG        mainRng;

    LightCacheSystem lightCacheSys_direct;
    LightCacheSystem lightCacheSys_indirect;

    struct InternalVariables {
        int currentRow;
        /* screen space detail */
        bool isScreenSpaceDetailComputed;
        /* direct lighting */
        bool isDirectLightingComputed;
        int  resolution_directLighting_cur;
        bool isDirectLightingKdTreeBuilt;
        /* indirect lighting */
        bool isIndirectLightingComputed;
        int  resolution_indirectLighting_cur;
        bool isIndirectLightingKdTreeBuilt;
        bool isIndirectLightingVisualized;
        /* final render */
        bool isFinalRenderCompleted;
        int  finalRender_curPass;
    } states;

protected:
    void saveLightCachePreview(const char* file);
    bool computeDirectLighting(InteractiveRenderer* renderer);
    bool buildKdTree_directLighting(InteractiveRenderer* renderer);
    bool computeIndirectLighting(InteractiveRenderer* renderer);
    bool buildKdTree_indirectLighting(InteractiveRenderer* renderer);
    bool visualizeIndirectLighting(InteractiveRenderer* renderer);
    bool finalRender(InteractiveRenderer* renderer);

    bool checkIfNeedToRefineLighting(int x, int y, int size);
    void computeScreenSpaceDetail(const RenderBuffer& albedoBuffer,
        const RenderBuffer& normalBuffer, const RenderBuffer& depthBuffer);

    VEC3 getIndirectLightingIntensity(const VEC3& p, const VEC3& n);
    VEC3 getDirectLightingIntensity(const VEC3& p, const VEC3& n);

    int  findValidLightCachePosition(ray& r, ray_query& Q);
    int  findValidRenderPosition(ray& r, ray_query& Q);

public:
    virtual void createRenderTask(
        Scene* scene, Camera* camera, RenderSettings* settings,
        RenderTask_Albedo* albedoTask,
        RenderTask_Normal* normalTask,
        RenderTask_Depth* depthTask
    );
    virtual void interactiveSubmit(InteractiveRenderer* renderer);
    virtual RenderStatus interactiveRender(InteractiveRenderer* renderer);
    virtual void onRenderStart(InteractiveRenderer* renderer);
    virtual void onRenderFinish(InteractiveRenderer* renderer);
    virtual bool checkRenderTask(String& log);

};
