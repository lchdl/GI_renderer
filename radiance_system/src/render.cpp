#include "render.h"

bool g_debug = false;


VEC3 computeReflectionDir(VEC3 n, VEC3 d) {
    return d - n * 2 * dot(n, d);
}
VEC3 computeRefractionDir(VEC3 d, VEC3 n, REAL ior_in, REAL ior_out) {
    REAL ior_rel = ior_in / ior_out;
    REAL cosTheta = dot(n, d);
    if (cosTheta >= 0) {
        n = -n;
        cosTheta = -cosTheta;
        if (ior_rel * Sqrt(ONE - cosTheta * cosTheta) > ONE) {
            return computeReflectionDir(n, d);
        }
    }
    VEC3 V = -d;
    REAL cosine = dot(n, V);
    return n * (ior_rel * cosine - Sqrt(ONE - ior_rel *
        ior_rel * (ONE - cosine * cosine))) - V * ior_rel;
}
REAL clampComp(const REAL& v, const REAL& compMin, const REAL& compMax)
{
    REAL o = v;
    if (o < compMin) o = compMin;
    if (o > compMax) o = compMax;
    return o;
}
VEC3 clampComp(const VEC3& v, const REAL& compMin, const REAL& compMax)
{
    VEC3 o = v;
    if (o.x < compMin) o.x = compMin;
    if (o.x > compMax) o.x = compMax;
    if (o.y < compMin) o.y = compMin;
    if (o.y > compMax) o.y = compMax;
    if (o.z < compMin) o.z = compMin;
    if (o.z > compMax) o.z = compMax;
    return o;
}
VEC4 clampComp(const VEC4& v, const REAL& compMin, const REAL& compMax)
{
    VEC4 o = v;
    if (o.s < compMin) o.s = compMin;
    if (o.s > compMax) o.s = compMax;
    if (o.x < compMin) o.x = compMin;
    if (o.x > compMax) o.x = compMax;
    if (o.y < compMin) o.y = compMin;
    if (o.y > compMax) o.y = compMax;
    if (o.z < compMin) o.z = compMin;
    if (o.z > compMax) o.z = compMax;
    return o;
}
VEC3 ARGBToRGB(const VEC4& v)
{
    return VEC3(v.x, v.y, v.z);
}
int ipow(int base, int power)
{
    int res = 1;
    while (power > 0) {
        res *= base;
        power--;
    }
    return res;
}
REAL computeSmoothness(const VEC3& n1, const VEC3& n2)
{
    /* return a value between [0,1] to indicate smoothness between two */
    /* surface normals. 0: not smooth at all (180deg), 1: coplanar (smooth) */
    return dot(n1, n2) * HALF + HALF;
}
void computeNormalMappingTBN(const VEC3& p1, const VEC3& p2, const VEC3& p3, const VEC2& t1, const VEC2& t2, const VEC2& t3, VEC3& T, VEC3& B, VEC3& N)
{
    /* https://learnopengl.com/Advanced-Lighting/Normal-Mapping */
    /* compute local basis for a triangle p1-p2-p3 with corresponding */
    /* texture coordinate t1-t2-t3 */
    VEC2 D1 = t2 - t1;
    VEC2 D2 = t3 - t2;
    const REAL& dU1 = D1.x;
    const REAL& dV1 = D1.y;
    const REAL& dU2 = D2.x;
    const REAL& dV2 = D2.y;
    REAL det = ONE / (dU1 * dV2 - dU2 * dV1);
    VEC3 E1 = p2 - p1;
    VEC3 E2 = p3 - p2;
    T = normalize(VEC3(
        det * (dV2 * E1.x - dV1 * E2.x),
        det * (dV2 * E1.y - dV1 * E2.y),
        det * (dV2 * E1.z - dV1 * E2.z)));
    B = normalize(VEC3(
        det * (-dU2 * E1.x + dU1 * E2.x),
        det * (-dU2 * E1.y + dU1 * E2.y),
        det * (-dU2 * E1.z + dU1 * E2.z)));
    N = normalize(cross(T, B));
}
void computeNormalMappingTBN(const triangleEx& tri, VEC3& T, VEC3& B, VEC3& N)
{
    VEC2 t1 = VEC2(tri.t[0].x, tri.t[0].y);
    VEC2 t2 = VEC2(tri.t[1].x, tri.t[1].y);
    VEC2 t3 = VEC2(tri.t[2].x, tri.t[2].y);
    return computeNormalMappingTBN(
        tri.p[0], tri.p[1], tri.p[2],
        t1, t2, t3,
        T, B, N);
}
void normalMappingTransformToWorld(const VEC3& T, const VEC3& B, const VEC3& N, const VEC3& localNormal, VEC3& worldNormal)
{
    worldNormal = MAT3x3(
        T.x, B.x, N.x,
        T.y, B.y, N.y,
        T.z, B.z, N.z)* localNormal;
}

/* utility function */
VEC3 rgbToVec3(BYTE r, BYTE g, BYTE b) {
    REAL denom = ONE / REAL(255.0);
    return VEC3(denom * r, denom * g, denom * b);
}
VEC4 rgbToVec4(BYTE r, BYTE g, BYTE b)
{
    REAL denom = ONE / REAL(255.0);
    return VEC4(ONE, denom * r, denom * g, denom * b);
}
VEC4 argbToVec4(BYTE a, BYTE r, BYTE g, BYTE b)
{
    REAL denom = ONE / REAL(255.0);
    return VEC4(denom * a, denom * r, denom * g, denom * b);
}
void albedoToRGBA(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a)
{
    *r = byteClamp(in->x * REAL(255.0));
    *g = byteClamp(in->y * REAL(255.0));
    *b = byteClamp(in->z * REAL(255.0));
    *a = byteClamp(in->s * REAL(255.0));
}
void colorToRGBA(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a)
{
    return albedoToRGBA(in, r, g, b, a);
}
void normalToRGBA(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a)
{
    *r = byteClamp((in->x + ONE) * REAL(128.0));
    *g = byteClamp((in->y + ONE) * REAL(128.0));
    *b = byteClamp((in->z + ONE) * REAL(128.0));
    *a = byteClamp(in->s * REAL(255.0));
}
void depthToRGBA(VEC4* in, BYTE* r, BYTE* g, BYTE* b, BYTE* a) {
    REAL scale = REAL(0.02);
    VEC3 inNorm = ONE - VEC3(
        Arctan(in->x * scale) * TWO * PI_INV,
        Arctan(in->y * scale) * TWO * PI_INV,
        Arctan(in->z * scale) * TWO * PI_INV);
    *r = byteClamp(inNorm.x * REAL(255.0));
    *g = byteClamp(inNorm.y * REAL(255.0));
    *b = byteClamp(inNorm.z * REAL(255.0));
    *a = byteClamp(in->s * REAL(255.0));
}
BYTE byteClamp(const REAL& value)
{
    if (value < 0) return BYTE(0);
    else if (value > 255) return BYTE(255);
    else return BYTE(value);
}
void vflipBuffer(BYTE* p, int stride, int h)
{
    BYTE* temp = (BYTE*)malloc(stride);
    int i = 0, j = h - 1;
    while (i < j) {
        memcpy(temp, p + i * stride, stride);
        memcpy((void*)(p + i * stride), p + j * stride, stride);
        memcpy((void*)(p + j * stride), temp, stride);
        i++, j--;
    }
    free(temp);
}

Bitmap::Bitmap()
{
    w = h = 0;
    c = 0;
    data = NULL;
}
Bitmap::~Bitmap()
{
    unload();
}
Bitmap::Bitmap(const Bitmap & that)
{
    this->w = that.w;
    this->h = that.h;
    this->c = that.c;
    
    int bytes = this->w * this->h * this->c;
    if (bytes > 0)
        this->data = malloc(bytes);
    else
        this->data = NULL;

    if (this->data != NULL)
        memcpy(this->data, that.data, bytes);
}
Bitmap Bitmap::operator=(const Bitmap & that)
{
    if (this == &that)
        return (*this);

    unload();

    this->w = that.w;
    this->h = that.h;
    this->c = that.c;

    int bytes = this->w * this->h * this->c;
    if (bytes > 0)
        this->data = malloc(bytes);
    else
        this->data = NULL;

    if (this->data != NULL)
        memcpy(this->data, that.data, bytes);

    return (*this);
}
bool Bitmap::load(const char * file, bool vflip)
{
    unload();
    
    if (vflip) stbi_set_flip_vertically_on_load(1);
    else stbi_set_flip_vertically_on_load(0);

    unsigned char* raw = stbi_load(file, &this->w, &this->h, &this->c, 0);
    if (raw == NULL) {
        printf("STB_IMAGE error: %s\n", stbi_failure_reason());
        return false;
    }
    if (this->c < 3 || this->c > 4) {
        printf("image with %d channels is unsupported.\n", this->c);
        stbi_image_free(raw);
        return false;
    }
    /* force 4 channel */
    this->data = (unsigned char*)malloc(this->w * this->h * 4);
    /* stbi load image as RGBA, we need to convert it to ARGB */
    if (this->c == 3) {
        /* RGB -> ARGB */
        unsigned char pixel[4];
        pixel[0] = 255; /* A=255 */
        unsigned char* srcptr = raw, 
            *dstptr = (unsigned char*)this->data;
        for (int y = 0; y < this->h; y++) {
            for (int x = 0; x < this->w; x++) {
                pixel[1] = srcptr[0];
                pixel[2] = srcptr[1];
                pixel[3] = srcptr[2];
                memcpy(dstptr, pixel, 4);
                dstptr += 4;
                srcptr += 3;
            }
        }
    }
    else if (this->c == 4) {
        /* RGBA -> ARGB */
        unsigned char pixel[4];
        unsigned char* srcptr = raw,
            *dstptr = (unsigned char*)this->data;
        for (int y = 0; y < this->h; y++) {
            for (int x = 0; x < this->w; x++) {
                pixel[0] = srcptr[3];
                pixel[1] = srcptr[0];
                pixel[2] = srcptr[1];
                pixel[3] = srcptr[2];
                /*if (pixel[0] > 0)
                    __debugbreak();*/
                memcpy(dstptr, pixel, 4);
                dstptr += 4;
                srcptr += 4;
            }
        }
    }
    /* image loaded */
    stbi_image_free(raw);
    return true;
}
bool Bitmap::save(const char * file, bool vflip)
{
    int bytes = 4 * w * h;
    unsigned char * pnew = (unsigned char*)malloc(bytes);
    int j = 0;
    unsigned char c;
    for (int i = 0; i < w*h; i++) {
        c = ((unsigned char*)data)[i * 4 + 1]; /* r */
        pnew[j++] = c;
        c = ((unsigned char*)data)[i * 4 + 2]; /* g */
        pnew[j++] = c;
        c = ((unsigned char*)data)[i * 4 + 3]; /* b */
        pnew[j++] = c;
        c = ((unsigned char*)data)[i * 4]; /* a */
        pnew[j++] = c;
    }
    if (vflip)
        vflipBuffer(pnew, w * 4, h);
    stbi_write_png(file, w, h, 4, pnew, w * 4);
    free(pnew);
    return true;
}
void Bitmap::create(int w, int h)
{
    unload();

    this->w = w;
    this->h = h;
    this->c = 4;
    this->data = malloc(w * h * this->c);

}
void Bitmap::unload()
{
    if (data) {
        free(data);
        data = NULL;
    }
    w = h = c = 0;
}
int Bitmap::getWidth() const
{
    return w;
}
int Bitmap::getHeight() const
{
    return h;
}
int Bitmap::toOffset(int x, int y)
{
    return y * this->w + x;
}
VEC4 Bitmap::toVec4(unsigned int x)
{
    /* little-endian */
    REAL _z = REAL((x & 0xff000000) >> 24);
    REAL _y = REAL((x & 0x00ff0000) >> 16);
    REAL _x = REAL((x & 0x0000ff00) >> 8);
    REAL _s = REAL(x & 0x000000ff);
    return VEC4(_s, _x, _y, _z);
}
VEC4 Bitmap::sample(VEC3 uvw)
{
/*

    x3---------x4
    |          |
dv  |          |
^   |          |
|   |          |
|   x1---------x2
|
+----> du

*/
    REAL u = wrap(uvw.x,ONE), v = wrap(uvw.y,ONE); /* ignore w component and wrap to [0,1) */

    REAL U = u * REAL(this->w), V = v * REAL(this->h);
    int X0 = int(U), Y0 = int(V);
    int X1 = X0 + 1, Y1 = Y0 + 1;
    REAL du = U - X0 , dv = V - Y0;

    if (X0 < 0) X0 = 0;
    if (X0 > this->w - 1) X0 = this->w - 1;
    if (Y0 < 0) Y0 = 0;
    if (Y0 > this->h - 1) Y0 = this->h - 1;

    if (X1 < 0) X1 = 0;
    if (X1 > this->w - 1) X1 = this->w - 1;
    if (Y1 < 0) Y1 = 0;
    if (Y1 > this->h - 1) Y1 = this->h - 1;

    unsigned int _x1 = ((unsigned int*)data)[toOffset(X0, Y0)];
    unsigned int _x2 = ((unsigned int*)data)[toOffset(X1, Y0)];
    unsigned int _x3 = ((unsigned int*)data)[toOffset(X0, Y1)];
    unsigned int _x4 = ((unsigned int*)data)[toOffset(X1, Y1)];

    REAL norm = ONE / REAL(255.0);

    VEC4 x1 = toVec4(_x1) * norm;
    VEC4 x2 = toVec4(_x2) * norm;
    VEC4 x3 = toVec4(_x3) * norm;
    VEC4 x4 = toVec4(_x4) * norm;

    VEC4 lo = x1 + (x2 - x1) * du;
    VEC4 hi = x3 + (x4 - x3) * du;
    VEC4 mi = lo + (hi - lo) * dv;

    return mi;
}
const char * Bitmap::getClassName() const
{
    return "Bitmap";
}
Bitmap Bitmap::stitchBitmaps(const Array<String>& files)
{
    if (files.size() == 0)
        return Bitmap();

    Bitmap result;

    int w = 0, h = 0;
    int bmx = 0, bmy = 0;
    int length = int(Sqrt(REAL(files.size())) + REAL(0.5));
    if (length * length < files.size()) length++;

    for (int i = 0; i < files.size(); i++) {
        Bitmap bm;
        if (bm.load(files[i].c_str(), false) == false) {
            return Bitmap();
        }
        if (w == 0 && h == 0) {
            w = bm.getWidth();
            h = bm.getHeight();
            result.create(w * length, h * length);
        }
        else {
            if (w != bm.getWidth() || h != bm.getHeight()) {
                printf("Error, bitmap size changed when stitching.\n");
                return Bitmap();
            }
        }
        for (int srcY = 0; srcY < h; srcY++) {
            for (int srcX = 0; srcX < w; srcX++) {
                int dstX = bmx * w + srcX, 
                    dstY = bmy * h + srcY;
                unsigned int* dstptr = (unsigned int*)(result.data);
                unsigned int* srcptr = (unsigned int*)(bm.data);
                dstptr[dstY * result.w + dstX] = srcptr[srcY * w + srcX];
            }
        }

        bmx++;
        if (bmx >= length) { bmy++; bmx = 0; }
    }

    return result;
}

void Entity::updateWorldBoundingBox()
{
    aabb localBoundingBox = getLocalBoundingBox();
    bboxWorld = this->xform.transformBoundingBoxToWorld(localBoundingBox);
}
Entity::Entity()
{
    materialId = 0;
    visible = true;
}
Entity::~Entity()
{
}
void Entity::useMaterial(const int& materialId)
{
    this->materialId = materialId;
}
int Entity::getMaterialId() const
{
    return materialId;
}
aabb Entity::getWorldBoundingBox()
{
    return this->bboxWorld;
}
bool Entity::isVisible() const
{
    return visible;
}
void Entity::setVisible(bool state)
{
    visible = state;
}
VEC3 Entity::getLocalLightSample(XorwowRNG& rng) const
{
    if (lightSamples.size() > 0) {
        int rnd = rng.uniform_int(0, lightSamples.size());
        return lightSamples[rnd];
    }
    else {
        return VEC3::zeros();
    }
}
void Entity::setTransform(const Transform& xform)
{
    this->xform = xform;
    updateWorldBoundingBox();
}
void Entity::setTransform(const REAL& scale, const QUAT& rotation, const VEC3& translation)
{
    this->xform.setTransform(scale, rotation, translation);
    updateWorldBoundingBox();
}
Transform Entity::getTransform() const
{
    return xform;
}

Array<ray_hit> MeshEntity::rayQuery(const ray* r, const REAL& searchDist)
{
    Array<ray_hit> hits;
    /* transform world ray to local ray */
    ray localRay;
    localRay.o = xform.transformPointToLocal(r->o);
    localRay.d = xform.transformVectorToLocal(r->d);
    this->doRaytracing(localRay, hits, searchDist / xform.getScale());
    /* transform ray trace info from local to world */
    for (int i = 0; i < hits.size(); i++) {
        hits[i].dist *= xform.getScale();
        hits[i].normal = xform.transformVectorToWorld(hits[i].normal);
        hits[i].point = xform.transformPointToWorld(hits[i].point);
    }
    return hits;
}
Array<ray_hit> MeshEntity::rayQuery(const ray * r, const Transform * proxyXform, 
    const REAL & searchDist)
{
    Array<ray_hit> hits;
    /* transform world ray to local ray, but using proxy's transform */
    ray localRay;
    localRay.o = proxyXform->transformPointToLocal(r->o);
    localRay.d = proxyXform->transformVectorToLocal(r->d);
    this->doRaytracing(localRay, hits, searchDist / proxyXform->getScale());
    /* transform ray trace info from local to world */
    for (int i = 0; i < hits.size(); i++) {
        hits[i].dist *= proxyXform->getScale();
        hits[i].normal = proxyXform->transformVectorToWorld(hits[i].normal);
        hits[i].point = proxyXform->transformPointToWorld(hits[i].point);
    }
    return hits;
}
bool MeshEntity::loadAsset(const char * file)
{
    bool success = BSTMesh::loadAsset(file);
    if (success)
        updateWorldBoundingBox();
    return success;
}
const char * MeshEntity::getClassName() const
{
    return "MeshEntity";
}
void MeshEntity::buildLightSamples(const int & count)
{
    /* we dont need to build light samples multiple times */
    if (lightSamples.size() == 0)
        lightSamples = this->sampleMesh(count);
}
VEC3 MeshEntity::getWorldLightSample(XorwowRNG & rng) const
{
    VEC3 localSample = getLocalLightSample(rng);
    return this->xform.transformPointToWorld(localSample);
}
VEC3 MeshEntity::getWorldLightSample(XorwowRNG & rng, const Transform * proxyXform) const
{
    VEC3 localSample = this->getLocalLightSample(rng);
    return proxyXform->transformPointToWorld(localSample);
}
VEC3 MeshEntity::normalQuery(const int& triangleId, const VEC3& uvw, const Scene* scene) const
{
    VEC3 nLocal = scene->getMaterial(this->materialId)->getLocalNormal(uvw, scene);
    VEC3 T, B, N;
    triangleEx tri = getTriangleEx(triangleId);
    /* we also need to rotate the triangle if this entity has set user transforms */
    tri.p[0] = xform.transformPointToWorld(tri.p[0]);
    tri.p[1] = xform.transformPointToWorld(tri.p[1]);
    tri.p[2] = xform.transformPointToWorld(tri.p[2]);
    /* compute tangent, bitangent, normal */
    computeNormalMappingTBN(tri, T, B, N);
    VEC3 nWorld;
    /* transform local space to world space */
    normalMappingTransformToWorld(T, B, N, nLocal, nWorld);
    return nWorld;
}
VEC3 MeshEntity::normalQuery(const int& proxyMaterialId, const Transform* proxyXform, 
    const int& triangleId, const VEC3& uvw, const Scene* scene) const
{
    VEC3 nLocal = scene->getMaterial(proxyMaterialId)->getLocalNormal(uvw, scene);
    VEC3 T, B, N;
    triangleEx tri = getTriangleEx(triangleId);
    /* we also need to rotate the triangle if this entity has set user transforms */
    tri.p[0] = proxyXform->transformPointToWorld(tri.p[0]);
    tri.p[1] = proxyXform->transformPointToWorld(tri.p[1]);
    tri.p[2] = proxyXform->transformPointToWorld(tri.p[2]);
    /* compute tangent, bitangent, normal */
    computeNormalMappingTBN(tri, T, B, N);
    VEC3 nWorld;
    /* transform local space to world space */
    normalMappingTransformToWorld(T, B, N, nLocal, nWorld);
    return nWorld;
}
aabb MeshEntity::getLocalBoundingBox() const
{
    if (this->bstRoot == NULL)
        return aabb();
    else return this->bstRoot->box;
}
MeshEntity::MeshEntity()
{
    materialId = 0;
    setTransform(ONE, QUAT::identity(), VEC3::zeros());
}
MeshEntity::~MeshEntity()
{
}

Scene::SceneBuildSettings Scene::defaultBuildSettings()
{
    SceneBuildSettings settings;
    settings.maxDepth = 16;
    settings.sahResolution = 16;
    settings.splitThreshold = 4;
    return settings;
}
void Scene::buildBVH_splitBox(const aabb * box, const char axis, const REAL & value, aabb * left, aabb * right)
{
    switch (axis) {
    case 'x':
        left->bmin = box->bmin;
        left->bmax = VEC3(value, box->bmax.y, box->bmax.z);
        right->bmin = VEC3(value, box->bmin.y, box->bmin.z);
        right->bmax = box->bmax;
        break;
    case 'y':
        left->bmin = box->bmin;
        left->bmax = VEC3(box->bmax.x, value, box->bmax.z);
        right->bmin = VEC3(box->bmin.x, value, box->bmin.z);
        right->bmax = box->bmax;
        break;
    case 'z':
        left->bmin = box->bmin;
        left->bmax = VEC3(box->bmax.x, box->bmax.y, value);
        right->bmin = VEC3(box->bmin.x, box->bmin.y, value);
        right->bmax = box->bmax;
        break;
    default:
        break;
    }
}
void Scene::buildBVH_sahCost(const Array<int>& idxs, const aabb & left, const aabb & right, const char & axis, const REAL & value, REAL * cost)
{
    int leftCount = 0, rightCount = 0;
    for (int i = 0; i < idxs.size(); i++) {
        aabb box = resources.entities[idxs[i]]->getWorldBoundingBox();
        if (isIntersect(box, left)) leftCount++;
        if (isIntersect(box, right)) rightCount++;
    }
    *cost = surfaceArea(left) * REAL(leftCount) + surfaceArea(right) * REAL(rightCount);
}
void Scene::buildBVH_recur(BVHNode * node, const Array<int>& idxs)
{
    bool shouldSplit = true;
    aabb lBox, rBox;
    Array<int> lIdxs, rIdxs;

    /* determine if this node should split */
    if (node->depth < this->bvhSettings.maxDepth &&
        idxs.size() > this->bvhSettings.splitThreshold)
    {
        /* find best split */
        VEC3 dxyz = node->box.bmax - node->box.bmin;

        /* init axes and splits (three directions) */
        int& sahResolution = this->bvhSettings.sahResolution;
        REAL* costs = new REAL[sahResolution * 3];
        REAL* splits = new REAL[sahResolution * 3];
        char* axes = new char[sahResolution * 3];
        for (int i = 0; i < sahResolution * 3; i++) {
            int D = i / sahResolution, S = i % sahResolution;
            REAL step = dxyz.e[D] / REAL(sahResolution + 1);
            axes[i] = 'x' + D;
            splits[i] = node->box.bmin.e[D] + REAL(S + 1) * step;
        }

        /* find minimum cost for all dimensions */
        for (int i = 0; i < sahResolution * 3; i++) {
            buildBVH_splitBox(&(node->box), axes[i], splits[i], &lBox, &rBox);
            buildBVH_sahCost(idxs, lBox, rBox, axes[i], splits[i], &(costs[i]));
        }
        int ax = argmin(costs, 3 * sahResolution);
        char axis = axes[ax];
        REAL split = splits[ax];


        buildBVH_splitBox(&(node->box), axis, split, &lBox, &rBox);
        for (int i = 0; i < idxs.size(); i++) {
            aabb box = resources.entities[idxs[i]]->getWorldBoundingBox();
            if (isIntersect(box, lBox)) lIdxs.append(idxs[i]);
            if (isIntersect(box, rBox)) rIdxs.append(idxs[i]);
        }

        /* calculate volume overlap */
        const REAL maxIncreaseRatio = REAL(0.40);
        REAL increaseRatio = ONE - REAL(lIdxs.size() + rIdxs.size()) / REAL(idxs.size());
        if (node->depth > this->bvhSettings.maxDepth * REAL(0.75) &&
            increaseRatio > maxIncreaseRatio) {
            shouldSplit = false;
        }
        else
        {
            shouldSplit = true;
        }

        delete[] costs;
        delete[] splits;
        delete[] axes;
    }
    else {
        shouldSplit = false;
    }

    /* * * split topNode into two child nodes or turn it to leaf topNode * * */

    if (shouldSplit) {
        /* build nodes */
        node->left = new BVHNode();
        node->right = new BVHNode();
        node->left->depth = node->depth + 1;
        node->right->depth = node->depth + 1;
        node->left->box = lBox;
        node->right->box = rBox;

        /* resursive build */
        buildBVH_recur(node->left, lIdxs);
        buildBVH_recur(node->right, rIdxs);
    }
    else {
        node->n = idxs.size();
        node->entityIds = (int*)malloc(sizeof(int)*(node->n));
        for (int i = 0; i < node->n; i++) {
            node->entityIds[i] = idxs[i];
        }
    }
}
void Scene::buildBVH_boxFit(const Array<int>& idxs, VEC3 * bmin, VEC3 * bmax)
{
    VEC3 _bmin = VEC3(REAL(REAL_MAX), REAL(REAL_MAX), REAL(REAL_MAX));
    VEC3 _bmax = VEC3(REAL(-REAL_MAX), REAL(-REAL_MAX), REAL(-REAL_MAX));

    for (int i = 0; i < idxs.size(); i++) {
        aabb box = resources.entities[idxs[i]]->getWorldBoundingBox();
        if (_bmin.x > box.bmin.x) _bmin.x = box.bmin.x;
        if (_bmin.y > box.bmin.y) _bmin.y = box.bmin.y;
        if (_bmin.z > box.bmin.z) _bmin.z = box.bmin.z;
        if (_bmax.x < box.bmax.x) _bmax.x = box.bmax.x;
        if (_bmax.y < box.bmax.y) _bmax.y = box.bmax.y;
        if (_bmax.z < box.bmax.z) _bmax.z = box.bmax.z;
    }

    *bmin = _bmin;
    *bmax = _bmax;
}
void Scene::buildBVH(const SceneBuildSettings & settings)
{
    clearBVH();
    
    Array<int> allEntityIds;
    for (int i = 0; i < getEntityCount(); i++)
        allEntityIds.append(i);

    bvhRoot = new BVHNode();
    bvhRoot->depth = 1;
    buildBVH_boxFit(allEntityIds, &(bvhRoot->box.bmin), &(bvhRoot->box.bmax));
    buildBVH_recur(bvhRoot, allEntityIds);
}
void Scene::clearBVH_recur(BVHNode * node)
{
    if (node == NULL) return;
    if (node->left) {
        clearBVH_recur(node->left);
        node->left = NULL;
    }
    if (node->right) {
        clearBVH_recur(node->right);
        node->right = NULL;
    }
    if (node->entityIds) {
        free(node->entityIds);
        node->entityIds = NULL;
    }
    delete node;
}
void Scene::clearBVH()
{
    clearBVH_recur(bvhRoot);
    bvhRoot = NULL;
}
bool Scene::checkVisibilityBetween(const VEC3& p1, const VEC3 p2)
{
    REAL dist = len(p1 - p2);
    ray r(p1, p2 - p1);
    ray_query Q = rayQuery(&r, dist);
    if (Q.numHits() == 0)
        return true;
    else
        return false;
}
bool Scene::checkVisibilityBetween(const VEC3& p1, const VEC3 p2, 
    const VEC3& n1, const VEC3& n2)
{
    REAL smoothness = computeSmoothness(n1, n2);
    if (smoothness < REAL(0.6)) {
        return checkVisibilityBetween(p1, p2);
    }
    else {
        REAL x = len(p2 - p1);
        REAL l = x / dot(n1, n2);
        REAL d = Sqrt(l * l - x * x);
        const REAL minOffset = REAL(0.05);
        if (d < minOffset) d = minOffset;
        return checkVisibilityBetween(p1 + n1 * d, p2 + n2 * d);
    }
}
int Scene::addMeshEntity(const char * file)
{
    MeshEntity* entity = new MeshEntity();
    if (entity->loadAsset(file) == false) {
        delete entity;
        printf("cannot load mesh asset \"%s\".\n", file);
        return -1;
    }
    else {
        resources.entities.append(entity);
        return resources.entities.size() - 1;
    }
}
int Scene::loadBitmap(const char * file)
{
    Bitmap* bitmap = new Bitmap();
    if (bitmap->load(file) == false) {
        delete bitmap;
        printf("cannot load bitmap \"%s\".\n", file);
        return -1;
    }
    else {
        resources.maps.append(bitmap);
        return resources.maps.size() - 1;
    }
}
int Scene::addMaterial(const Material* material)
{
    Material* mat = new Material();
    *mat = *material;
    resources.materials.append(mat);
    return resources.materials.size() - 1;
}
int Scene::addCameraPreset(const Camera * camera)
{
    cameraPresets.append(*camera);
    return cameraPresets.size() - 1;
}
bool Scene::setEntityMaterial(const int & entityId, const int & materialId)
{
    if (entityId < 0 || entityId >= getEntityCount())
        return false;
    if (materialId < 0 || materialId >= getMaterialCount())
        return false;
    resources.entities[entityId]->useMaterial(materialId);
    Material* material = getMaterial(materialId);
    if (material->getMaterialType() & Material_Emissive) {
        getEntity(entityId)->buildLightSamples(material->getNumEmissiveSamples());
    }
    return true;
}
void Scene::setEntityTransform(const int & entityId, const Transform & xform)
{
    getEntity(entityId)->setTransform(xform);
}
void Scene::setEntityTransform(const int & entityId, const VEC3& translation, const QUAT& rotation, const REAL & scale)
{
    getEntity(entityId)->setTransform(scale, rotation, translation);
}
int Scene::copyEntity(const int & entityId)
{
    EntityProxy* proxy = new EntityProxy();
    Entity* proxyTarget = this->getEntity(entityId);
    proxy->setProxyTarget(proxyTarget);
    resources.entities.append(proxy);
    return resources.entities.size() - 1;
}
int Scene::getEntityCount() const
{
    return resources.entities.size();
}
const Camera& Scene::getCameraPreset(const int & cameraId) const
{
    return cameraPresets[cameraId];
}
Entity * Scene::getEntity(const int& entityId) const
{
    return resources.entities[entityId];
}
Material * Scene::getEntityMaterial(const int & entityId) const
{
    return resources.materials[getEntity(entityId)->getMaterialId()];
}
int Scene::getMaterialCount() const
{
    return resources.materials.size();
}
Material * Scene::getMaterial(const int& materialId) const
{
    return resources.materials[materialId];
}
Array<int> Scene::getEmissiveEntityIds() const
{
    return resources.emissiveEntityIds;
}
const SceneResources & Scene::getResources() const
{
    return resources;
}
bool Scene::buildScene(const SceneBuildSettings& settings)
{
    this->bvhSettings = settings;
    this->buildBVH(this->bvhSettings);
    /* check how many entities are actually light (with emissive material) */
    resources.emissiveEntityIds.clear();
    for (int i = 0; i < getEntityCount(); i++) {
        Material* material = getEntityMaterial(i);
        if (material->getMaterialType() & Material_Emissive) {
            /* make this entity light */
            resources.emissiveEntityIds.append(i);
            getEntity(i)->buildLightSamples(material->getNumEmissiveSamples());
        }
    }

    return true;
}
void Scene::unload()
{
    clearBVH();

    for (int i = 0; i < resources.entities.size(); i++)
        delete resources.entities[i];
    resources.entities.clear();

    for (int i = 0; i < resources.maps.size(); i++)
        delete resources.maps[i];
    resources.maps.clear();

    for (int i = 0; i < resources.materials.size(); i++)
        delete resources.materials[i];
    resources.materials.clear();
}
ray_query Scene::rayQuery(const ray* r, const REAL searchDist) const
{
    /* shoot single ray into the scene */
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BVHNode*> nodeStack;
    nodeStack.reserve(bvhSettings.maxDepth + 1);
    nodeStack.push(bvhRoot); /* initialize stack and trace */

    ray_query Q;
    int hit = 0; /* if anything in the scene is hit by this ray */

    while (nodeStack.isEmpty() == false) {
        BVHNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                int entityId = topNode->entityIds[i];
                Entity* entity = resources.entities[entityId];
                Material* material = getEntityMaterial(entityId);
                aabb box = entity->getWorldBoundingBox();
                if (isIntersect(*r, box)) {
                    Array<ray_hit> hits = entity->rayQuery(r, searchDist);
                    for (int i = 0; i < hits.size(); i++) {
                        hits[i].entityId = entityId;
                        Q.addHit(hits[i]);
                    }
                }
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(*r, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(*r, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
    Q.sort();
    return Q;
}
ray_query Scene::rayQueryNormalMapping(const ray* r, const REAL searchDist) const
{
    /* shoot single ray into the scene */
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BVHNode*> nodeStack;
    nodeStack.reserve(bvhSettings.maxDepth + 1);
    nodeStack.push(bvhRoot); /* initialize stack and trace */

    ray_query Q;
    int hit = 0; /* if anything in the scene is hit by this ray */

    while (nodeStack.isEmpty() == false) {
        BVHNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                int entityId = topNode->entityIds[i];
                Entity* entity = resources.entities[entityId];
                Material* material = getEntityMaterial(entityId);
                aabb box = entity->getWorldBoundingBox();
                if (isIntersect(*r, box)) {
                    Array<ray_hit> hits = entity->rayQuery(r, searchDist);
                    for (int i = 0; i < hits.size(); i++) {
                        hits[i].entityId = entityId;
                        /* if normal map is used, perform a second query */
                        if (material->getNormalMapId() >= 0)
                            hits[i].normal = entity->normalQuery(hits[i].triangleId, hits[i].uvw, this);
                        Q.addHit(hits[i]);
                    }
                }
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(*r, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(*r, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
    Q.sort();
    return Q;
}
ray_query Scene::rayQueryIgnoreMaterial(const ray * r, const unsigned int & ignoredMaterial, const REAL searchDist) const
{
    /* shoot single ray into the scene */
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BVHNode*> nodeStack;
    nodeStack.reserve(bvhSettings.maxDepth + 1);
    nodeStack.push(bvhRoot); /* initialize stack and trace */

    ray_query Q;
    int hit = 0; /* if anything in the scene is hit by this ray */

    while (nodeStack.isEmpty() == false) {
        BVHNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                int entityId = topNode->entityIds[i];
                Entity* entity = resources.entities[entityId];
                Material* material = getEntityMaterial(entityId);
                if (material->getMaterialType() & ignoredMaterial)
                    continue;
                aabb box = entity->getWorldBoundingBox();
                if (isIntersect(*r, box)) {
                    Array<ray_hit> hits = entity->rayQuery(r, searchDist);
                    for (int i = 0; i < hits.size(); i++) {
                        hits[i].entityId = entityId;
                        Q.addHit(hits[i]);
                    }
                }
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(*r, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(*r, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
    Q.sort();
    return Q;
}
ray_query Scene::rayQueryIgnoreEntity(const ray * r, const unsigned int & ignoredEntityId, const REAL searchDist) const
{
    /* shoot single ray into the scene */
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BVHNode*> nodeStack;
    nodeStack.reserve(bvhSettings.maxDepth + 1);
    nodeStack.push(bvhRoot); /* initialize stack and trace */

    ray_query Q;
    int hit = 0; /* if anything in the scene is hit by this ray */

    while (nodeStack.isEmpty() == false) {
        BVHNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                int entityId = topNode->entityIds[i];
                if (entityId == ignoredEntityId)
                    continue;
                Entity* entity = getEntity(entityId);
                Material* material = getEntityMaterial(entityId);
                aabb box = entity->getWorldBoundingBox();
                if (isIntersect(*r, box)) {
                    Array<ray_hit> hits = entity->rayQuery(r, searchDist);
                    for (int i = 0; i < hits.size(); i++) {
                        hits[i].entityId = entityId;
                        Q.addHit(hits[i]);
                    }
                }
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(*r, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(*r, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
    Q.sort();
    return Q;
}
VEC3 Scene::computeDirectLightIntensity(
    const VEC3 & samplePoint, const VEC3& normal, const int & samplesPerLight, XorwowRNG& rng)
{
    VEC3 averageIntensity = VEC3::zeros();
    int numLights = this->resources.emissiveEntityIds.size();
    for (int i = 0; i < numLights; i++) {
        int lightEntityId = this->resources.emissiveEntityIds[i];
        Entity* entity = getEntity(lightEntityId);
        for (int sampleId = 0; sampleId < samplesPerLight; sampleId++) {
            VEC3 localLightSample = entity->getLocalLightSample(rng);
            Transform& entityWorldTransform = entity->getTransform();
            VEC3 worldLightSample = entityWorldTransform.transformPointToWorld(localLightSample);
            VEC3 rayPath = worldLightSample - samplePoint;
            REAL d = len(rayPath); /* distance between light and sample point */
            ray r(samplePoint, rayPath);
            ray_query Q = rayQuery(&r);

            VEC3 accumulateThinFalloff = VEC3::ones();
            for (int i = 0; i < Q.numHits(); i++) {
                Material* material = getEntityMaterial(Q[i].entityId);
                if (material->getMaterialType() & Material_Diffuse) {
                    if (material->isThinObject())
                        accumulateThinFalloff *= material->getThinFalloff();
                    else
                        accumulateThinFalloff = VEC3::zeros();
                }
                else if (Q[i].entityId == lightEntityId) {
                    /* ray reached this light */
                    /* sample light color */
                    VEC4 emissiveIntensity = material->getEmissiveIntensity(Q[i].uvw, this);
                    /* compute falloff using smoothed inverse square root */
                    REAL falloff = material->getEmissiveFalloff();
                    REAL attenuation = ONE / (falloff * d * d + ONE); /* 1 / (fd^2 + 1) */
                    VEC3 finalSampleIntensity = Fabs(dot(r.d, normal)) * attenuation * emissiveIntensity.xyz() * accumulateThinFalloff;
                    averageIntensity += finalSampleIntensity;
                    break; /* quit loop once we hit the first light */
                }
                else {
                    /* hit non diffusive objects */
                    accumulateThinFalloff = VEC3::zeros();
                }
            }

        }
    }
    averageIntensity /= REAL(numLights * samplesPerLight);
    return averageIntensity;
}
VEC3 Scene::finalGather(const VEC3& samplePoint, const VEC3& normal, 
    const int& hemiSamples, const int& k, const LightCacheSystem& lightCacheSys, 
    const bool& checkVisibility, XorwowRNG& rng)
{
    VEC3 indirectIntensity = VEC3::zeros();

    for (int i = 0; i < hemiSamples; i++)
    {
        ray r;
        r.o = samplePoint;
        r.d = randHemi(normal, rng);
        ray_query Q = rayQuery(&r);

        if (Q.numHits() > 0) {
            if (dot(Q[0].normal, r.d) > ZERO)
                Q[0].normal = -Q[0].normal;
            VEC3 q = Q[0].point + Q[0].normal * REAL(0.001);
            Array<light_cache> C = lightCacheSys.knnSearch(q, k, false);
            VEC3 sampAvgInten = VEC3::zeros();
            int acceptedCache = 0;
            for (int i = 0; i < C.size(); i++) {
                if (checkVisibilityBetween(q, C[i].p) == true
                    && dot(Q[0].normal, C[i].n) > REAL(0.7)) {
                    /* only collect coplanar sample points */
                    sampAvgInten += C[i].i;
                    acceptedCache++;
                }
            }
            if (acceptedCache > 0)
                sampAvgInten /= REAL(acceptedCache);

            indirectIntensity += sampAvgInten;
        }
    }
    indirectIntensity /= REAL(hemiSamples);

    return indirectIntensity;





    //VEC3 indirectIntensity = VEC3::zeros();
    //int acceptedSample = 0;
    //for (int i = 0; i < hemiSamples; i++)
    //{
    //    ray r;
    //    r.o = samplePoint;
    //    r.d = randHemi(normal, rng);
    //    ray_query Q = rayQuery(&r);

    //    if (Q.numHits() > 0) {
    //        VEC3 q = Q[0].point;
    //        if (dot(Q[0].normal, r.d) > ZERO)
    //            Q[0].normal = -Q[0].normal;
    //        q = q + Q[0].normal * REAL(0.001);

    //        if (g_debug) {
    //            printf("\n");
    //            printVec3("      position:", q);
    //        }

    //        Array<light_cache> C = lightCacheSys.knnSearch(q, k, false);
    //        for (int i = 0; i < C.size(); i++) {
    //            if (checkVisibilityBetween(q, C[i].p) == true
    //                && dot(Q[0].normal, C[i].n) > REAL(0.7)) {
    //                /* only collect coplanar sample points */
    //                indirectIntensity += C[i].i;
    //                acceptedSample++;

    //                if (g_debug) {
    //                    printf("         accept:(%.3f,%.3f,%.3f), inten(%.3f,%.3f,%.3f)\n",
    //                        C[i].p.x, C[i].p.y, C[i].p.z, C[i].i.x, C[i].i.y, C[i].i.z);
    //                }

    //            }
    //            else {
    //                if (g_debug) {
    //                    printVec3("         reject:", C[i].p);
    //                }
    //            }
    //        }
    //    }
    //}
    //indirectIntensity /= REAL(acceptedSample);

    //return indirectIntensity;

}
Scene::Scene()
{
    bvhRoot = NULL;
}
Scene::~Scene()
{
    unload();
}

Material::Material()
{
    setDiffuse(VEC3(REAL(0.8), REAL(0.8), REAL(0.8)), -1, ONE, -1);
}
Material::~Material()
{
}
unsigned int Material::getMaterialType() const
{
    if (lensq(emissiveIntensity) > ZERO)
        return Material_Emissive;
    unsigned int mask = 0;
    if (lensq(diffuseColor) > ZERO)
        mask |= Material_Diffuse;
    if (lensq(reflectColor) > ZERO)
        mask |= Material_Reflect;
    if (lensq(refractColor) > ZERO)
        mask |= Material_Refract;
    if (lensq(volumeColor) > ZERO)
        mask |= Material_Volumetric;
    return mask;
}
VEC4 Material::getDiffuseColor(VEC3 uvw, const Scene* scene) const
{
    VEC4 texColor;
    if (diffuseMapId < 0)
        texColor = VEC4::ones();
    else
        texColor = scene->getResources().maps[diffuseMapId]->sample(uvw);
    return texColor * VEC4(ONE, diffuseColor);
}
VEC3 Material::getLocalNormal(VEC3 uvw, const Scene* scene) const
{
    if (normalMapId < 0)
        return VEC3(0, 0, 1);
    VEC4 normalColor = scene->getResources().maps[normalMapId]->sample(uvw);
    VEC3 n;
    n.x = TWO * (normalColor.x - HALF);
    n.y = TWO * (normalColor.y - HALF);
    n.z = TWO * (normalColor.z - HALF);
    return n;
}
int Material::getNormalMapId() const
{
    return normalMapId;
}
VEC4 Material::getEmissiveIntensity(VEC3 uvw, const Scene * scene) const
{
    VEC4 texColor;
    if (emissiveMapId < 0)
        texColor = VEC4::ones();
    else
        texColor = scene->getResources().maps[emissiveMapId]->sample(uvw);
    return texColor * VEC4(ONE, emissiveIntensity);
}
REAL Material::getEmissiveFalloff() const
{
    return emissiveFalloff;
}
int Material::getNumEmissiveSamples() const
{
    return emissiveSamplesNum;
}
REAL Material::getRefractionIndex() const
{
    return refractIndex;
}
void Material::setEmissive(VEC3 emissiveColor, int emissiveMapId, REAL emissiveFalloff, int emissiveSamples)
{
    this->emissiveIntensity = clampComp(emissiveColor, ZERO, REAL_MAX);
    this->emissiveMapId = emissiveMapId;
    this->emissiveFalloff = clampComp(emissiveFalloff, ZERO, REAL_MAX);
    this->emissiveSamplesNum = emissiveSamples;
}
void Material::setDiffuse(VEC3 diffuseColor, int diffuseMapId, REAL diffuseRoughness, 
    int normalMapId, bool thinObject, VEC3 thinObjectFalloff)
{
    this->diffuseColor = clampComp(diffuseColor, ZERO, ONE);
    this->diffuseMapId = diffuseMapId;
    this->diffuseRoughness = clampComp(diffuseRoughness, ZERO, ONE);
    this->normalMapId = normalMapId;
    this->thinObject = thinObject;
    this->thinObjectFalloff = thinObjectFalloff;
}
void Material::setReflect(VEC3 reflectColor, int reflectMapId, REAL reflectRoughness, bool fresnelReflect)
{
    this->reflectColor = clampComp(reflectColor, ZERO, ONE);
    this->reflectMapId = reflectMapId;
    this->reflectRoughness = clampComp(reflectRoughness, ZERO, ONE);
    this->fresnelReflect = fresnelReflect;
}
void Material::setRefract(VEC3 refractColor, REAL refractIndex, REAL refractRoughness)
{
    this->refractColor = clampComp(refractColor, ZERO, ONE);
    this->refractIndex = clampComp(refractIndex, ZERO, REAL_MAX);
    this->refractRoughness = clampComp(refractRoughness, ZERO, ONE);
}
void Material::setVolumetric(VEC3 volumeColor, REAL volumeFalloff)
{
    this->volumeColor = clampComp(volumeColor, ZERO, ONE);
    this->volumeFalloff = clampComp(volumeFalloff, ZERO, REAL_MAX);
}
bool Material::hasProperty(unsigned int materialTypeMask)
{
    unsigned int type = getMaterialType();
    if (type & materialTypeMask) return true;
    else return false;
}

bool Material::isThinObject() const
{
    return this->thinObject;
}

VEC3 Material::getThinFalloff() const
{
    return this->thinObjectFalloff;
}

Scene::BVHNode::BVHNode()
{
    depth = 0;
    n = 0;
    entityIds = NULL;
    left = right = NULL;
}
Scene::BVHNode::~BVHNode()
{
    if (left) { delete left; left = NULL; }
    if (right) { delete right; right = NULL; }
    if (entityIds) { free(entityIds); entityIds = NULL; }
}
bool Scene::BVHNode::isLeaf() const
{
    return (left == NULL || right == NULL);
}
Scene::SceneBuildSettings::SceneBuildSettings()
{
    maxDepth = 16;
    splitThreshold = 4;
    sahResolution = 16;
}
bool Scene::SceneBuildSettings::isValid()
{
    if (maxDepth < 1) return false;
    if (splitThreshold < 1) return false;
    if (sahResolution < 1) return false;
    return true;
}

RenderBuffer::RenderBuffer()
{
    w = h = 0;
    buffer = NULL;
    raw = NULL;
}
RenderBuffer::~RenderBuffer()
{
    destroy();
}
RenderBuffer::RenderBuffer(const RenderBuffer & that)
{
    this->w = that.w;
    this->h = that.h;
    int bufferSize = sizeof(VEC4) * w * h;
    this->buffer = (VEC4*)malloc(bufferSize);
    this->raw = (BYTE*)malloc(sizeof(BYTE) * 4 * w * h);
    memcpy(this->buffer, that.buffer, bufferSize);
    memcpy(this->raw, that.raw, sizeof(BYTE) * 4 * w * h);
}
void RenderBuffer::destroy()
{
    if (buffer) {
        free(buffer);
        buffer = NULL;
    }
    if (raw) {
        free(raw);
        raw = NULL;
    }
    w = h = 0;
}
void RenderBuffer::create(int w, int h)
{
    destroy();

    this->w = w;
    this->h = h;

    int bufferSize = sizeof(VEC4) * w * h;
    buffer = (VEC4*)malloc(bufferSize);
    memset(buffer, 0, bufferSize);
    this->raw = (BYTE*)malloc(sizeof(BYTE) * 4 * w * h);
}
void RenderBuffer::setPixel(int x, int y, const VEC4 & value)
{
    if (x < 0 || x >= this->w || y < 0 || y >= this->h)
        return;
    int offset = y * this->w + x;
    buffer[offset] = value;
}
void RenderBuffer::setBlock(int x, int y, int w, int h, const VEC4 & value)
{
    for (int y0 = y; y0 < y + h; y0++) {
        for (int x0 = x; x0 < x + w; x0++) {
            if (x0 < 0 || x0 >= getWidth() || y0 < 0 || y0 >= getHeight())
                continue;
            buffer[y0 * getWidth() + x0] = value;
        }
    }
}
VEC4 RenderBuffer::getPixel(int x, int y) const
{
    if (x < 0 || x >= this->w || y < 0 || y >= this->h)
        return VEC4();
    int offset = y * this->w + x;
    return buffer[offset];
}
int RenderBuffer::getWidth() const
{
    return this->w;
}
int RenderBuffer::getHeight() const
{
    return this->h;
}
RenderBuffer RenderBuffer::operator+(const RenderBuffer & that)
{
    RenderBuffer rb;
    rb.create(w, h);
    for (int i = 0; i < w*h; i++)
        rb.buffer[i] = this->buffer[i] + that.buffer[i];
    return rb;
}
void RenderBuffer::operator+=(const RenderBuffer & that)
{
    for (int i = 0; i < w*h; i++)
        buffer[i] = buffer[i] + that.buffer[i];
}
RenderBuffer RenderBuffer::operator*(const REAL & v)
{
    RenderBuffer rb;
    rb.create(w, h);
    for (int i = 0; i < w*h; i++)
        rb.buffer[i] = this->buffer[i] * v;
    return rb;
}
RenderBuffer RenderBuffer::operator/(const REAL & v)
{
    REAL v_inv = ONE / v;
    RenderBuffer rb;
    rb.create(w, h);
    for (int i = 0; i < w*h; i++)
        rb.buffer[i] = this->buffer[i] * v_inv;
    return rb;
}
void RenderBuffer::operator/=(const REAL & v)
{
    REAL v_inv = ONE / v;
    for (int i = 0; i < w*h; i++)
        buffer[i] = buffer[i] * v_inv;
}
RenderBuffer RenderBuffer::operator=(const RenderBuffer & that)
{
    if (this == &that)
        return (*this);
    
    destroy();

    this->w = that.w;
    this->h = that.h;
    int bufferSize = sizeof(VEC4) * w * h;
    this->buffer = (VEC4*)malloc(bufferSize);
    memcpy(this->buffer, that.buffer, bufferSize);
    this->raw = (BYTE*)malloc(sizeof(BYTE) * 4 * w * h);
    memcpy(this->raw, that.raw, sizeof(BYTE) * 4 * w * h);

    return (*this);
}
void RenderBuffer::saveBitmap(pfnColorFunction func, const char * file)
{
    /* save image to disk (RGBA 8 bits per channel) */
    stbi_write_png(file, w, h, 4, asBitmap(func), w * 4);
}
void RenderBuffer::saveBitmap(pfnColorFunction func, REAL scale, const char * file)
{
    /* save image to disk (RGBA 8 bits per channel) */
    stbi_write_png(file, w, h, 4, asBitmap(func, scale), w * 4);
}
const BYTE * RenderBuffer::asBitmap(pfnColorFunction func)
{
    /* save image to disk (RGBA 8 bits per channel) */
    for (int i = 0; i < w * h; i++) {
        BYTE r, g, b, a;
        func(buffer + i, &r, &g, &b, &a);
        raw[i * 4] = r;
        raw[i * 4 + 1] = g;
        raw[i * 4 + 2] = b;
        raw[i * 4 + 3] = a;
    }
    return raw;
}
const BYTE * RenderBuffer::asBitmap(pfnColorFunction func, REAL scale)
{
    for (int i = 0; i < w * h; i++) {
        BYTE r, g, b, a;
        VEC4 px = buffer[i];
        px.x *= scale;
        px.y *= scale;
        px.z *= scale;
        func(&px, &r, &g, &b, &a);
        raw[i * 4] = r;
        raw[i * 4 + 1] = g;
        raw[i * 4 + 2] = b;
        raw[i * 4 + 3] = a;
    }
    return raw;
}
void RenderBuffer::collectRegion(int x, int y, int w, int h, Array<VEC4>& pixels, bool keepShape) const
{
    if (keepShape == false) {
        for (int y0 = y; y0 < y + h; y0++) {
            for (int x0 = x; x0 < x + w; x0++) {
                if (x0 < 0 || x0 >= this->w || y0 < 0 || y0 >= this->h)
                    continue;
                pixels.append(buffer[y0 * this->w + x0]);
            }
        }
    }
    else {
        for (int y0 = y; y0 < y + h; y0++) {
            for (int x0 = x; x0 < x + w; x0++) {
                if (x0 < 0 || x0 >= this->w || y0 < 0 || y0 >= this->h)
                    pixels.append(VEC4());
                else 
                    pixels.append(buffer[y0 * this->w + x0]);
            }
        }
    }
}
VEC4 RenderBuffer::max() const
{
    VEC4 m = VEC4::min();
    for (int i = 0; i < w * h; i++) {
        VEC4 px = this->buffer[i];
        if (m.x < px.x) m.x = px.x;
        if (m.y < px.y) m.y = px.y;
        if (m.z < px.z) m.z = px.z;
    }
    return m;
}

void InteractiveRenderer::logMessage(const char * msg, const int& level)
{
    LogItem item;
    if (logs.isFull())
        logs.get(item);
    item.msg = msg;
    item.level = level;
    logs.put(item);
}
void InteractiveRenderer::initFrontendDisplay(int w, int h)
{
    const char* fontPath = "res/sys.ttf";
    const int fontSize = 16;
    const int logBufferSize = 64;
    mainViewport.gen(0, 0, w, h);
    mainViewport.setWindowMargin(w, h);
    sysFont.loadFont(fontPath, fontSize);
    sysFont.setWindowMargin(w, h);
    logs.create(logBufferSize);
    texBlit.gen();
    texBlit.setWindowMargin(w, h);
    this->windowWidth = w;
    this->windowHeight = h;

    logInfo("** GI rendering system, developed by Chenghao Liu **");
}
void InteractiveRenderer::onWindowResize(int w, int h)
{
    /* save "real" window size (w and h could be zero) */
    this->windowWidth = w;
    this->windowHeight = h;
    /* force w and h to be at least 1 */
    if (w <= 0) w = 1;
    if (h <= 0) h = 1;
    mainViewport.setViewportMargin(0.0f, 0.0f, float(w), float(h));
    mainViewport.setWindowMargin(w, h);
    sysFont.setWindowMargin(w, h);
    texBlit.setWindowMargin(w, h);
}
void InteractiveRenderer::logInfo(const char * msg, bool suppressDuplicates)
{
    int level = 0;
    if (logs.sizeUsed() > 0) {
        int last = logs.sizeUsed() - 1;
        if (logs[last].msg == msg && suppressDuplicates)
            return;
        else
            logMessage(msg, level);
    }
    else
        logMessage(msg, level);
}
void InteractiveRenderer::logWarning(const char * msg, bool suppressDuplicates)
{
    int level = 1;
    if (logs.sizeUsed() > 0) {
        int last = logs.sizeUsed() - 1;
        if (logs[last].msg == msg && suppressDuplicates)
            return;
        else
            logMessage(msg, level);
    }
    else
        logMessage(msg, level);
}
void InteractiveRenderer::logError(const char * msg, bool suppressDuplicates)
{
    int level = 2;
    if (logs.sizeUsed() > 0) {
        int last = logs.sizeUsed() - 1;
        if (logs[last].msg == msg && suppressDuplicates)
            return;
        else
            logMessage(msg, level);
    }
    else
        logMessage(msg, level);
}
void InteractiveRenderer::updateLog(const char * msg, const int& level)
{
    if (logs.sizeUsed() == 0) {
        LogItem item;
        item.level = level;
        item.msg = msg;
        logs.put(item);
    }
    else {
        int last = logs.sizeUsed() - 1;
        logs[last].msg = msg;
        logs[last].level = level;
    }
}
int InteractiveRenderer::getWindowWidth() const
{
    return this->windowWidth;
}
int InteractiveRenderer::getWindowHeight() const
{
    return this->windowHeight;
}
InteractiveRenderer::InteractiveRenderer()
{
    prevStatus = Render_NoTask;
    renderCalls = 0;
    windowWidth = windowHeight = 0;
    fatalError = false;
}
InteractiveRenderer::~InteractiveRenderer()
{
}
void InteractiveRenderer::display()
{
    mainViewport.bind();
    mainViewport.clear();

    glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    
    /* display canvas */
    if (getCurrentRenderTask()!=NULL) {
        if (getCurrentRenderTask()->displayTexture.isInit() == false) {
            logError("ERROR, display texture for the render task is not ready! "
                "You need to manually initialize \"displayTexture\" when render task is created.");
        }
        int w = getCurrentRenderTask()->displayTexture.getWidth();
        int h = getCurrentRenderTask()->displayTexture.getHeight();
        int x = (getWindowWidth() - w) / 2;
        int y = (getWindowHeight() - h) / 2;
        texBlit.draw(getCurrentRenderTask()->displayTexture, float(x), float(y));
    }

    /* display logs */
    REAL cursorX = REAL(3.0), cursorY = REAL(3.0);
    REAL textScale = REAL(1.0);
    for (int i = logs.sizeUsed() - 1; i >= 0; i--) {
        VEC3 color;
        switch (logs[i].level) {
        case 0: /* INFO */
            color = VEC3::all(REAL(0.9));
            break;
        case 1:
            color = VEC3(ONE, ONE, ZERO);
            break;
        case 2:
            color = VEC3(ONE, ZERO, ZERO);
            break;
        default:
            color = VEC3(ONE, ZERO, ZERO);
            break;
        }
        sysFont.drawText(logs[i].msg.c_str(), cursorX, cursorY, textScale, color);
        cursorY += REAL(10.0);
    }

    mainViewport.unbind();
    mainViewport.submit();
}
void InteractiveRenderer::addRenderTask(RenderTask * task)
{
    if (renderTasks.sizeMax() == 0)
        renderTasks.create(64);
    if (renderTasks.put(task) == false)
        logError("ERROR, render queue is full.");
}
bool InteractiveRenderer::checkRenderTasks()
{
    if (this->renderTasks.sizeUsed() == 0) {
        logInfo("INFO: no render task.");
        return true;
    }
    bool ok = true;
    for (int i = 0; i < this->renderTasks.sizeUsed(); i++) {
        RenderTask* renderTask = this->renderTasks[i];
        String log;
        if (renderTask->checkRenderTask(log) == false) {
            char buf[64];
            sprintf(buf, "WARNING in task [%d]: ", i);
            String s = buf;
            s += log;
            logWarning(s.c_str());
            ok = false;
        }
    }
    return ok;
}
RenderTask * InteractiveRenderer::getCurrentRenderTask() const
{
    if (renderTasks.sizeUsed() == 0)
        return NULL;
    else
        return renderTasks[0];
}
RenderStatus InteractiveRenderer::interactiveRender() {
    
    if (fatalError) {
        display();
        return Render_NoTask;
    }

    RenderStatus currentStatus;

    if (getCurrentRenderTask() != NULL) {
        if (renderCalls == 0)
            getCurrentRenderTask()->onRenderStart(this);
        currentStatus = getCurrentRenderTask()->interactiveRender(this);
        renderCalls++;
        if (prevStatus != Render_Finished && currentStatus == Render_Finished)
            getCurrentRenderTask()->onRenderFinish(this);
        getCurrentRenderTask()->interactiveSubmit(this);

        if (prevStatus != Render_Finished && currentStatus == Render_Finished) {
            RenderTask* task;
            renderTasks.get(task);
            prevStatus = Render_NoTask;
            renderCalls = 0;
        }
        else {
            prevStatus = currentStatus;
        }
    }
    else {
        currentStatus = Render_NoTask;
        logInfo("INFO: render queue is empty.");
    }
    display();
    return currentStatus;
}
void InteractiveRenderer::setErrorFlagAndQuit(const String& errMsg)
{
    fatalError = true;
    logError(errMsg.c_str());
}

int ray_query::numHits() const
{
    return hits.size();
}
ray_hit& ray_query::operator[](const int & idx)
{
    return hits[idx];
}
void ray_query::clear()
{
    hits.clear();
}
void ray_query::removeFirstN(const int & N)
{
    ray_query Qnew;
    for (int i = N; i < numHits(); i++) {
        Qnew.addHit((*this)[i]);
    }
    hits = Qnew.hits;
}
void ray_query::sort()
{
    int numRecords = hits.size();
    for (int i = 0; i < numRecords; i++) {
        for (int j = i+1; j < numRecords; j++) {
            if (hits[i].dist > hits[j].dist) {
                swap(hits[i], hits[j]);
            }
        }
    }
}
void ray_query::addHit(ray_hit& hit)
{
    hits.append(hit);
}

unsigned int RenderTask::toRGBA(BYTE r, BYTE g, BYTE b)
{
    return ((((unsigned int)(r)) << 24) | (((unsigned int)(g)) << 16) | (((unsigned int)(b)) << 8) | 0xFF);
}
unsigned int RenderTask::toRGBA(BYTE r, BYTE g, BYTE b, BYTE a)
{
    return ((((unsigned int)(r)) << 24) | (((unsigned int)(g)) << 16) | (((unsigned int)(b)) << 8) | a);
}
void RenderTask::drawRect(BYTE * p, int x, int y, int w, int h, unsigned int color, int W, int H)
{
    for (int _x = x; _x < x + w; _x++)
        drawPoint(p, _x, y, color, W, H);
    for (int _x = x; _x < x + w; _x++)
        drawPoint(p, _x, y + h - 1, color, W, H);
    for (int _y = y; _y < y + h; _y++)
        drawPoint(p, x, _y, color, W, H);
    for (int _y = y; _y < y + h; _y++)
        drawPoint(p, x + w - 1, _y, color, W, H);
}
void RenderTask::drawPoint(BYTE * p, int x, int y, unsigned int color, int W, int H)
{
    unsigned int* pi = (unsigned int*)p;
    if (x < 0 || x >= W || y < 0 || y >= H)
        return;
    pi[y*W + x] = color;
}
String RenderTask::formatTime(double sec)
{
    char buf[64];
    if (sec < 0) {
        sprintf(buf, "N/A");
    }
    else if (sec < 1) {
        sprintf(buf, "%.3lfs", sec);
    }
    else if (sec < 10) {
        sprintf(buf, "%.2lfs", sec);
    }
    else if (sec < 30) {
        sprintf(buf, "%.1lfs", sec);
    }
    else if (sec < 60) {
        sprintf(buf, "%ds", int(sec));
    }
    else if (sec < 60 * 60) {
        sprintf(buf, "%dm,%ds", int(sec) / 60, int(sec) % 60);
    }
    else if (sec < 24 * 60 * 60) {
        sprintf(buf, "%dh,%dm", int(sec) / 3600, (int(sec) % 3600) / 60);
    }
    else if (sec < 30 * 24 * 60 * 60) {
        sprintf(buf, "%dd,%dh", int(sec) / (3600 * 24), (int(sec) % (3600*24)) / 3600);
    }
    else {
        sprintf(buf, ">30d");
    }
    return String(buf);
}
String RenderTask::formatProgress(const int & total, const int & finished)
{
    char buf[256];
    double time_elapse = timer.query();
    double time_remain = time_elapse * REAL(total - finished) / REAL(finished);
    String rtime = formatTime(time_remain);
    String etime = formatTime(time_elapse);
    sprintf(buf, "Progress: %.2lf%% | Elapsed: %s | Remaining: %s",
        double(finished) / double(total) * 100.0,
        etime.c_str(), rtime.c_str());
    return buf;
}
String RenderTask::formatSize(int bytes)
{
    char buf[64];
    if (bytes < 0) {
        sprintf(buf, "N/A");
    }
    else if (bytes == 1) {
        sprintf(buf, "1byte");
    }
    else if (bytes < 1024) {
        sprintf(buf, "%dbytes", bytes);
    }
    else if (bytes < 1024 * 1024) {
        sprintf(buf, "%.2fkB", REAL(bytes) / REAL(1024));
    }
    else if (bytes < 1024 * 1024 * 1024) {
        sprintf(buf, "%.2fMB", REAL(bytes) / REAL(1024 * 1024));
    }
    else {
        sprintf(buf, "%.2fGB", REAL(bytes) / REAL(1024 * 1024 * 1024));
    }
    return String(buf);
}
void RenderTask::onRenderStart(InteractiveRenderer * renderer) {
    timer.tick();
}
void RenderTask::onRenderFinish(InteractiveRenderer * renderer) 
{
}
bool RenderTask::checkRenderTask(String & log)
{
    return true;
}

EntityStack::EntityStack()
{
    this->reserve(8);
}
EntityStack::~EntityStack()
{
}
void EntityStack::pushEntityId(const int & entityId)
{
    this->push(entityId);
}
int EntityStack::popEntityId()
{
    if (this->isEmpty()) return -1;
    else {
        int i;
        this->pop(i);
        return i;
    }
}
REAL EntityStack::getCurrentIOR(const Scene* scene) const
{
    if (this->isEmpty()) return ONE; /* air */
    else {
        return scene->getEntityMaterial(top())->getRefractionIndex();
    }
}

void RenderTask_Albedo::saveImage(const char * file) {
    renderBuffer.saveBitmap(albedoToRGBA, file);
}
void RenderTask_Albedo::createRenderTask(Scene * scene, Camera * camera, RenderSettings * settings, const char * save)
{
    this->scene = scene;
    this->camera = *camera;
    this->settings = *settings;
    this->savePath = save;
    states.scan_y = 0;
    renderBuffer.create(settings->render.w, settings->render.h);
    displayTexture.gen(settings->render.w, settings->render.h, GL_RGBA, GL_NEAREST);
}
RenderStatus RenderTask_Albedo::interactiveRender(InteractiveRenderer * renderer) 
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;

    /* render one row at a time */
    const int rows = getNumProcessors();
    int y_start = states.scan_y;
    int y_end = states.scan_y + rows > h ? h : states.scan_y + rows;

    /* start rendering for each pixel */
#pragma omp parallel for schedule(static)
    for (int y = y_start; y < y_end; y++) {
        for (int x = 0; x < w; x++) {
            VEC4 pixelColor;
            ray r = camera.raycast(w, h, x, y);
            while (true) {
                ray_query Q = scene->rayQuery(&r);
                if (Q.numHits() > 0) {
                    Material* material = scene->getEntityMaterial(Q[0].entityId);
                    pixelColor = material->getDiffuseColor(Q[0].uvw, scene);
                    if (pixelColor.s < HALF) /* object texture has opacity */
                    {
                        if (dot(r.d, Q[0].normal) < ZERO)
                            Q[0].normal = -Q[0].normal;
                        r.o = Q[0].point + Q[0].normal * REAL(0.001);                        
                    }
                    else
                        break;
                }
                else {
                    pixelColor = VEC4::zeros();
                    break;
                }
            }
            renderBuffer.setPixel(x, y, pixelColor);
        }
    }
    /* update scanline, so next time it will continue */
    states.scan_y = y_end;

    /* update message */
    String progress = formatProgress(h, y_end);
    renderer->updateLog(progress.c_str());

    if (y_end >= h)
        return RenderStatus::Render_Finished;
    else
        return RenderStatus::Render_Running;
}
void RenderTask_Albedo::interactiveSubmit(InteractiveRenderer * renderer)
{
    /* submit rendered result to GPU texture, so that it can be seen */
    /* from the users. */
    const BYTE* p = this->renderBuffer.asBitmap(albedoToRGBA);
    int stride = sizeof(BYTE) * 4 * this->renderBuffer.getWidth();
    int w = this->renderBuffer.getWidth();
    int h = this->renderBuffer.getHeight();
    int bytes = stride * h;
    BYTE* pc = (BYTE*)malloc(bytes);
    memcpy(pc, p, bytes);
    /* draw annotations */
    if (states.scan_y < h) memset(pc + states.scan_y * stride, 255, stride);
    /* GL texture is vertically flipped, so we need to flip them too */
    vflipBuffer(pc, stride, h);
    displayTexture.update(pc);
    free(pc);
}
void RenderTask_Albedo::onRenderStart(InteractiveRenderer * renderer)
{
    RenderTask::onRenderStart(renderer);

    renderer->logInfo("Start rendering albedo map...");
    renderer->logInfo("");
}
void RenderTask_Albedo::onRenderFinish(InteractiveRenderer * renderer)
{
    RenderTask::onRenderFinish(renderer);

    renderer->logInfo("Render finished.");
    double dt = timer.tick();
    char buffer[512];
    sprintf(buffer, "Render time: %.2lfsec.", dt);
    renderer->logInfo(buffer);
    if (savePath.length() > 0)
    {
        saveImage(savePath.c_str());
        sprintf(buffer, "Rendered image saved to \"%s\".", savePath.c_str());
        renderer->logInfo(buffer);
    }
}
RenderBuffer RenderTask_Albedo::getRenderBuffer() const
{
    return renderBuffer;
}

Transform::Transform()
{
    scale = ONE;
    rotation = QUAT::identity();
    translation = VEC3::zeros();
}
Transform::Transform(const REAL & scale, const QUAT & rotation, const VEC3 & translation)
{
    setTransform(scale, rotation, translation);
}
Transform::~Transform()
{
}
Transform Transform::operator=(const Transform & that)
{
    this->scale = that.scale;
    this->rotation = normalize(that.rotation);
    this->translation = that.translation;
    return (*this);
}
void Transform::setTransform(const REAL& scale, const QUAT & rotation, const VEC3 & translation)
{
    this->scale = scale;
    this->rotation = normalize(rotation);
    this->translation = translation;
}
void Transform::setTransform(const Transform & xform)
{
    this->operator=(xform);
}
VEC3 Transform::transformPointToLocal(const VEC3 & p) const
{
    /* transform order: scale <- rotate <- translate */
    VEC3 pLocal = p;
    pLocal -= translation;
    pLocal = quatToMatrix(inv(rotation)) * pLocal;
    pLocal /= scale;
    return pLocal;
}
Array<VEC3> Transform::transformPointToLocal(const Array<VEC3>& p) const
{
    Array<VEC3> P;
    MAT3x3 rotationMatrix = quatToMatrix(inv(rotation));
    for (int i = 0; i < p.size(); i++) {
        VEC3 pLocal = p[i];
        pLocal -= translation;
        pLocal = rotationMatrix * pLocal;
        pLocal /= scale;
        P.append(pLocal);
    }
    return P;
}
VEC3 Transform::transformVectorToLocal(const VEC3 & v) const
{
    VEC3 vLocal;
    vLocal = quatToMatrix(inv(rotation)) * v;
    return vLocal;
}
Array<VEC3> Transform::transformVectorToLocal(const Array<VEC3>& v) const
{
    Array<VEC3> V;
    MAT3x3 rotationMatrix = quatToMatrix(inv(rotation));
    for (int i = 0; i < v.size(); i++) {
        VEC3 vLocal = v[i];
        vLocal = rotationMatrix * vLocal;
        V.append(vLocal);
    }
    return V;
}
VEC3 Transform::transformPointToWorld(const VEC3 & p) const
{
    /* transform order: scale -> rotate -> translate */
    VEC3 pWorld = p;
    pWorld *= scale;
    pWorld = quatToMatrix(rotation) * pWorld;
    pWorld += translation;
    return pWorld;
}
Array<VEC3> Transform::transformPointToWorld(const Array<VEC3>& p) const
{
    Array<VEC3> P;
    MAT3x3 rotationMatrix = quatToMatrix(rotation);
    for (int i = 0; i < p.size(); i++) {
        VEC3 pLocal = p[i];
        pLocal *= scale;
        pLocal = rotationMatrix * pLocal;
        pLocal += translation;
        P.append(pLocal);
    }
    return P;
}
VEC3 Transform::transformVectorToWorld(const VEC3 & v) const
{
    VEC3 vWorld;
    vWorld = quatToMatrix(rotation) * v;
    return vWorld;
}
Array<VEC3> Transform::transformVectorToWorld(const Array<VEC3>& v) const
{
    Array<VEC3> V;
    MAT3x3 rotationMatrix = quatToMatrix(rotation);
    for (int i = 0; i < v.size(); i++) {
        VEC3 vWorld;
        vWorld = rotationMatrix * vWorld;
        V.append(vWorld);
    }
    return V;
}
aabb Transform::transformBoundingBoxToWorld(const aabb & bbox) const
{
    Array<VEC3> pointsLocal;
    pointsLocal.reserve(8);

    REAL xmin = bbox.bmin.x, ymin = bbox.bmin.y, zmin = bbox.bmin.z;
    REAL xmax = bbox.bmax.x, ymax = bbox.bmax.y, zmax = bbox.bmax.z;

    pointsLocal.append(VEC3(xmin, ymin, zmin));
    pointsLocal.append(VEC3(xmax, ymin, zmin));
    pointsLocal.append(VEC3(xmin, ymax, zmin));
    pointsLocal.append(VEC3(xmax, ymax, zmin));
    pointsLocal.append(VEC3(xmin, ymin, zmax));
    pointsLocal.append(VEC3(xmax, ymin, zmax));
    pointsLocal.append(VEC3(xmin, ymax, zmax));
    pointsLocal.append(VEC3(xmax, ymax, zmax));

    Array<VEC3> pointsWorld = transformPointToWorld(pointsLocal);

    return computeBoundingBoxFromPointCloud(pointsWorld);
}
REAL Transform::getScale() const
{
    return scale;
}
QUAT Transform::getRotation() const
{
    return rotation;
}
VEC3 Transform::getTranslation() const
{
    return translation;
}

EntityProxy::EntityProxy()
{
    proxyTarget = NULL;
}
EntityProxy::~EntityProxy()
{
}
void EntityProxy::setProxyTarget(Entity * entity)
{
    proxyTarget = entity;
    this->materialId = entity->getMaterialId();
}
VEC3 EntityProxy::normalQuery(const int& triangleId, const VEC3& uvw, const Scene* scene) const
{
    return proxyTarget->normalQuery(this->materialId, &(this->xform) ,triangleId, uvw, scene);
}
VEC3 EntityProxy::normalQuery(const int& proxyMaterialId, const Transform* proxyXform, const int& triangleId, const VEC3& uvw, const Scene* scene) const
{
    return proxyTarget->normalQuery(proxyMaterialId, proxyXform, triangleId, uvw, scene);
}
void EntityProxy::buildLightSamples(const int & count)
{
    this->proxyTarget->buildLightSamples(count);
}
VEC3 EntityProxy::getLocalLightSample(XorwowRNG & rng) const
{
    return this->proxyTarget->getLocalLightSample(rng);
}
VEC3 EntityProxy::getWorldLightSample(XorwowRNG & rng) const
{
    return this->proxyTarget->getWorldLightSample(rng, &this->xform);
}
VEC3 EntityProxy::getWorldLightSample(XorwowRNG & rng, const Transform * proxyXform) const
{
    return this->proxyTarget->getWorldLightSample(rng, proxyXform);
}
aabb EntityProxy::getLocalBoundingBox() const
{
    return this->proxyTarget->getLocalBoundingBox();
}
const char * EntityProxy::getClassName() const
{
    return "EntityProxy";
}
Array<ray_hit> EntityProxy::rayQuery(const ray * r, const REAL & searchDist)
{
    return this->proxyTarget->rayQuery(r, &this->xform, searchDist);
}
Array<ray_hit> EntityProxy::rayQuery(const ray * r, const Transform * proxyXform, const REAL & searchDist)
{
    return this->proxyTarget->rayQuery(r, proxyXform, searchDist);
}

void RenderTask_Normal::saveImage(const char * file)
{
    renderBuffer.saveBitmap(normalToRGBA, file);
}
RenderStatus RenderTask_Normal::interactiveRender(InteractiveRenderer * renderer)
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;

    /* render one row at a time */
    const int rows = getNumProcessors();
    int y_start = states.scan_y;
    int y_end = states.scan_y + rows > h ? h : states.scan_y + rows;

    /* start rendering for each pixel */
#pragma omp parallel for schedule(static)
    for (int y = y_start; y < y_end; y++) {
        for (int x = 0; x < w; x++) {
            VEC4 normalColor;
            ray r = camera.raycast(w, h, x, y);
            while (true) {
                ray_query Q = scene->rayQueryNormalMapping(&r);
                if (Q.numHits() > 0) {
                    Material* material = scene->getEntityMaterial(Q[0].entityId);
                    normalColor = VEC4(ONE, Q[0].normal);
                    VEC4 pixelColor = material->getDiffuseColor(Q[0].uvw, scene);
                    if (pixelColor.s < HALF) /* object texture has opacity */
                    {
                        if (dot(r.d, Q[0].normal) < ZERO) {
                            r.o = Q[0].point - Q[0].normal * REAL(0.001);
                        }
                        else {
                            r.o = Q[0].point + Q[0].normal * REAL(0.001);
                        }
                    }
                    else
                        break;
                }
                else{
                    normalColor = VEC4::zeros();
                    break;
                }
            }

            renderBuffer.setPixel(x, y, normalColor);
        }
    }
    /* update scanline, so next time it will continue */
    states.scan_y = y_end;

    /* update message */
    String progress = formatProgress(h, y_end);
    renderer->updateLog(progress.c_str());

    if (y_end >= h)
        return RenderStatus::Render_Finished;
    else
        return RenderStatus::Render_Running;
}
void RenderTask_Normal::interactiveSubmit(InteractiveRenderer * renderer)
{    
    /* submit rendered result to GPU texture, so that it can be seen */
    /* from the users. */
    const BYTE* p = this->renderBuffer.asBitmap(normalToRGBA);
    int stride = sizeof(BYTE) * 4 * this->renderBuffer.getWidth();
    int w = this->renderBuffer.getWidth();
    int h = this->renderBuffer.getHeight();
    int bytes = stride * h;
    BYTE* pc = (BYTE*)malloc(bytes);
    memcpy(pc, p, bytes);
    /* draw annotations */
    if (states.scan_y < h) memset(pc + states.scan_y * stride, 255, stride);
    /* GL texture is vertically flipped, so we need to flip them too */
    vflipBuffer(pc, stride, h);
    displayTexture.update(pc);
    free(pc);
}
void RenderTask_Normal::onRenderStart(InteractiveRenderer * renderer)
{
    RenderTask::onRenderStart(renderer);

    renderer->logInfo("Start rendering normal map...");
    renderer->logInfo("");
}

void RenderTask_Depth::saveImage(const char * file)
{
    renderBuffer.saveBitmap(depthToRGBA, file);
}
RenderStatus RenderTask_Depth::interactiveRender(InteractiveRenderer * renderer)
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;

    /* render one row at a time */
    const int rows = getNumProcessors();
    int y_start = states.scan_y;
    int y_end = states.scan_y + rows > h ? h : states.scan_y + rows;

    /* start rendering for each pixel */
#pragma omp parallel for schedule(static)
    for (int y = y_start; y < y_end; y++) {
        for (int x = 0; x < w; x++) {
            VEC4 pixelColor;
            ray r = camera.raycast(w, h, x, y);
            while (true) {
                ray_query Q = scene->rayQuery(&r);
                if (Q.numHits() > 0) {
                    REAL depth = len(Q[0].point - camera.pos);

                    Material* material = scene->getEntityMaterial(Q[0].entityId);
                    VEC4 albedo = material->getDiffuseColor(Q[0].uvw, scene);
                    if (albedo.s < HALF) /* object texture has opacity */
                    {
                        if (dot(r.d, Q[0].normal) < ZERO) {
                            r.o = Q[0].point - Q[0].normal * REAL(0.001);
                        }
                        else {
                            r.o = Q[0].point + Q[0].normal * REAL(0.001);
                        }
                    }
                    else {
                        pixelColor = VEC4(ONE, VEC3::all(depth));
                        break;
                    }
                }
                else {
                    pixelColor = VEC4::zeros();
                    break;
                }
            }
            renderBuffer.setPixel(x, y, pixelColor);
        }
    }
    /* update scanline, so next time it will continue */
    states.scan_y = y_end;

    /* update message */
    String progress = formatProgress(h, y_end);
    renderer->updateLog(progress.c_str());

    if (y_end >= h)
        return RenderStatus::Render_Finished;
    else
        return RenderStatus::Render_Running;
}
void RenderTask_Depth::interactiveSubmit(InteractiveRenderer * renderer)
{
    /* submit rendered result to GPU texture, so that it can be seen */
    /* from the users. */
    const BYTE* p = this->renderBuffer.asBitmap(depthToRGBA);
    int stride = sizeof(BYTE) * 4 * this->renderBuffer.getWidth();
    int w = this->renderBuffer.getWidth();
    int h = this->renderBuffer.getHeight();
    int bytes = stride * h;
    BYTE* pc = (BYTE*)malloc(bytes);
    memcpy(pc, p, bytes);
    /* draw annotations */
    if (states.scan_y < h) memset(pc + states.scan_y * stride, 255, stride);
    /* GL texture is vertically flipped, so we need to flip them too */
    vflipBuffer(pc, stride, h);
    displayTexture.update(pc);
    free(pc);
}
void RenderTask_Depth::onRenderStart(InteractiveRenderer * renderer)
{
    RenderTask::onRenderStart(renderer);

    renderer->logInfo("Start rendering depth map...");
    renderer->logInfo("");
}

void RenderSettings::validateParameters()
{
    if (render.w <= 1) render.w = 1;
    if (render.h <= 1) render.h = 1;
    if (GI.traceDepth <= 1) GI.traceDepth = 1;
    if (GI.resolution_directLighting_min > 0) GI.resolution_directLighting_min = 0;
    if (GI.numSamplesPerLight < 1) GI.numSamplesPerLight = 1;
    if (GI.knnSamples_directLighting < 1) GI.knnSamples_directLighting = 1;
    if (GI.knnSamples_indirectLighting < 1) GI.knnSamples_indirectLighting = 1;
    if (GI.maxKdTreeDepth < 1) GI.maxKdTreeDepth = 1;
    if (GI.numFinalGatherSamples < 1) GI.numFinalGatherSamples = 1;
    if (GI.resolution_indirectLighting_min > 0)
        GI.resolution_indirectLighting_min = 0;
    if (GI.resolution_indirectLighting_max > 0)
        GI.resolution_indirectLighting_max = 0;
    if (GI.resolution_indirectLighting_min > GI.resolution_indirectLighting_max)
        GI.resolution_indirectLighting_min = GI.resolution_indirectLighting_max;
    if (GI.resolution_directLighting_min > 0)
        GI.resolution_directLighting_min = 0;
    if (GI.resolution_directLighting_max > 0)
        GI.resolution_directLighting_max = 0;
    if (GI.resolution_directLighting_min > GI.resolution_directLighting_max)
        GI.resolution_directLighting_min = GI.resolution_directLighting_max;
}
RenderSettings RenderSettings::defaultRenderSettings()
{
    RenderSettings settings;
    settings.render.w = 640;
    settings.render.h = 480;
    settings.render.antialiasing = true;
    settings.render.backgroundColor = VEC4::zeros();
    settings.GI.randomSampling = false;
    settings.GI.traceDepth = 3;
    settings.GI.numSamplesPerLight = 128;
    settings.GI.knnSamples_directLighting = 40;
    settings.GI.knnSamples_indirectLighting = 50;
    settings.GI.maxKdTreeDepth = 24;
    settings.GI.numFinalGatherSamples = 200;
    settings.GI.resolution_directLighting_min = -2;
    settings.GI.resolution_directLighting_max = -1;
    settings.GI.resolution_indirectLighting_min = -3;
    settings.GI.resolution_indirectLighting_max = -1;
    settings.GI.checkVisibility = true;

    return settings;
}

void RenderTask_GI::saveLightCachePreview(const char * file)
{
    renderBuffer.saveBitmap(albedoToRGBA, file);
}
bool RenderTask_GI::computeDirectLighting(InteractiveRenderer* renderer)
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;
    const int& numSamplesPerLight = settings.GI.numSamplesPerLight;
    const int step = ipow(2, -states.resolution_directLighting_cur);

    int totalRows = (h + step - 1) / step;
    int row = getNumProcessors();
    int rowStart = states.currentRow;
    int rowEnd = rowStart + row > totalRows ? totalRows : rowStart + row;

    /* uniform scatter in screen space and bounce */
#pragma omp parallel for schedule(static)
    for (int rowId = rowStart; rowId < rowEnd; rowId++) {
        int y0 = rowId * step;
        XorwowRNG rng(y0, 166);
        Array<light_cache> rowLightCache;
        for (int x0 = 0; x0 < w; x0 += step) {

            REAL x = REAL(x0) + REAL(step) * HALF, y = REAL(y0) + REAL(step) * HALF;
            ray r = camera.raycast(w, h, x0, y0, settings.GI.randomSampling, step, &rng);
            int localSize = step * 8 - 1;
            if ((states.resolution_directLighting_cur != settings.GI.resolution_directLighting_min)
                && checkIfNeedToRefineLighting(int(x), int(y), localSize) == false)
            {
                continue;
            }

            Array<light_cache> pxLightCache;
            pxLightCache.reserve(settings.GI.traceDepth);
            int bounce = settings.GI.traceDepth;

            VEC3 visualizeIntensity;

            while (bounce > 0) {
                ray_query Q;
                int status = findValidLightCachePosition(r, Q);
                if (status == -1) 
                {
                    break;
                }
                else if (status == 0) 
                {
                    VEC3 n = Q[0].normal;
                    if (dot(n, r.d) > ZERO)
                        n = -n; /* reverse normal if surface is a thin plane */
                    VEC3 p = Q[0].point + n * REAL(0.001);
                    VEC3 d = randHemi(n, rng);
                    /* set next status and continue */
                    r.o = p;
                    r.d = d;

                    /* visualization */
                    Material* material = scene->getEntityMaterial(Q[0].entityId);
                    if (bounce == settings.GI.traceDepth &&
                        (material->getMaterialType() & Material_Emissive)) 
                    {
                        visualizeIntensity = material->getEmissiveIntensity(Q[0].uvw, scene).xyz();
                    }
                }
                else if (status == 1) 
                {
                    Material* material = scene->getEntityMaterial(Q[0].entityId);
                    VEC3 n = Q[0].normal;
                    if (dot(n, r.d) > ZERO)
                        n = -n; /* reverse normal if surface is a thin plane */
                    VEC3 p = Q[0].point + n * REAL(0.001);
                    VEC3 d = randHemi(n, rng);
                    r.o = p;
                    r.d = d;
                    light_cache cache;
                    cache.p = p;
                    cache.n = n;
                    cache.i = ARGBToRGB(material->getDiffuseColor(Q[0].uvw, scene));
                    cache.i *= scene->computeDirectLightIntensity(cache.p, cache.n, numSamplesPerLight, rng);
                    pxLightCache.append(cache);
                    if (bounce == settings.GI.traceDepth)
                        visualizeIntensity = cache.i;
                }
                bounce--;
            }

            rowLightCache.append(pxLightCache);
            renderBuffer.setBlock(x0, y0, step, step, VEC4(ONE, visualizeIntensity));

        }
#pragma omp critical(RenderTask_GI_computeLightCache)
        lightCacheSys_direct.addLightCache(rowLightCache);
    }
    states.currentRow = rowEnd;

    String progressString = formatProgress(totalRows, rowEnd);
    char buf[512] = { 0 };
    sprintf(buf, "Current [%d]: ", states.resolution_directLighting_cur);
    String message = buf;
    message += progressString;
    renderer->updateLog(message);

    if (rowEnd == totalRows)
    {
        states.currentRow = 0;
        states.resolution_directLighting_cur++;
        if (states.resolution_directLighting_cur > settings.GI.resolution_directLighting_max) {
            renderer->logInfo("Finished.");
            String s = "Light cache computation finished in ";
            s += formatTime(timer.tick());
            s += ".";
            renderer->logInfo(s.c_str());
            char buf[256];
            sprintf(buf, "Constructing kd tree (%d samples), this may take a while...", lightCacheSys_direct.numCaches());
            renderer->logInfo(buf);
            if (lightCacheSys_direct.numCaches() > 200000) {
                /* it takes 20 minute to build a light cache with 1,000,000 samples */
                renderer->logWarning("WARNING: too many samples (> 200,000), building process might be slow");
            }
            /* visualize light cache position */
            //for (int i = 0; i < lightCacheSys_direct.numCaches(); i++) {
            //    int sx, sy;
            //    camera.worldPositionToScreen(lightCacheSys_direct.getCache(i).p,
            //        settings.render.w, settings.render.h, sx, sy);
            //    renderBuffer.setPixel(sx, sy, VEC4(ONE, ONE, ZERO, ZERO));
            //}
            saveLightCachePreview("direct.png");
            return true;
        }
        else {
            renderer->logInfo("");
            timer.tick();
            return false;
        }
    }
    else
        return false;
}
bool RenderTask_GI::buildKdTree_directLighting(InteractiveRenderer* renderer)
{
    Timer timer;
    timer.tick();
    lightCacheSys_direct.buildKdTree(settings.GI.knnSamples_directLighting, settings.GI.maxKdTreeDepth);
    String dt = formatTime(timer.tick());
    String s = "kd-tree construction finished in ";
    s += dt;
    s += ".";
    renderer->logInfo(s.c_str());
    s = "kd-tree memory consumption: ";
    s += formatSize(lightCacheSys_direct.computeCacheSize());
    s += ".";
    renderer->logInfo(s.c_str());
    renderer->logInfo("");
    return true;
}
bool RenderTask_GI::computeIndirectLighting(InteractiveRenderer* renderer)
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;
    //const int& numSamplesPerLight = settings.GI.numSamplesPerLight;
    const int step = ipow(2, -states.resolution_indirectLighting_cur);

    int totalRows = (h + step - 1) / step;
    int row = getNumProcessors();
    int rowStart = states.currentRow;
    int rowEnd = rowStart + row > totalRows ? totalRows : rowStart + row;

    /* uniform scatter in screen space and bounce */
#pragma omp parallel for schedule(static)
    for (int rowId = rowStart; rowId < rowEnd; rowId++) {
        int y0 = rowId * step;
        XorwowRNG rng(y0, 166);
        Array<light_cache> rowLightCache;
        for (int x0 = 0; x0 < w; x0 += step) {
            if (x0 == 24 && y0 == 582) {
                g_debug = true;
            }
            else
            {
                g_debug = false;
            }
            REAL x = REAL(x0) + REAL(step) * HALF, y = REAL(y0) + REAL(step) * HALF;
            ray r = camera.raycast(w, h, x0, y0, settings.GI.randomSampling, step, &rng);

            int localSize = step * 8 - 1;
            if ((states.resolution_indirectLighting_cur != settings.GI.resolution_indirectLighting_min)
                && checkIfNeedToRefineLighting(int(x),int(y), localSize) == false)
            {
                continue;
            }

            Array<light_cache> pxLightCache;
            pxLightCache.reserve(settings.GI.traceDepth);
            int bounce = settings.GI.traceDepth;
            
            VEC3 visualizeIntensity;

            while (bounce > 0) {
                ray_query Q;
                int status = findValidLightCachePosition(r, Q);
                if (status == -1)
                {
                    break;
                }
                else if (status == 0)
                {
                    VEC3 n = Q[0].normal;
                    VEC3 p = Q[0].point + n * REAL(0.001);
                    if (dot(n, r.d) > ZERO)
                        n = -n; /* reverse normal if surface is a thin plane */
                    VEC3 d = randHemi(n, rng);
                    /* set next status and continue */
                    r.o = p;
                    r.d = d;
                    /* visualization */
                    Material* material = scene->getEntityMaterial(Q[0].entityId);
                    if (bounce == settings.GI.traceDepth &&
                        (material->getMaterialType() & Material_Emissive))
                    {
                        visualizeIntensity = material->getEmissiveIntensity(Q[0].uvw, scene).xyz();
                    }
                }
                else if (status == 1)
                {
                    VEC3 n = Q[0].normal;
                    if (dot(n, r.d) > ZERO)
                        n = -n;
                    VEC3 p = Q[0].point + n * REAL(0.001);
                    VEC3 d = randHemi(n, rng);
                    r.o = p;
                    r.d = d;

                    light_cache cache;
                    cache.i = scene->finalGather(p, n,
                        settings.GI.numFinalGatherSamples,
                        settings.GI.knnSamples_directLighting,
                        lightCacheSys_direct, settings.GI.checkVisibility,
                        rng);
                    cache.p = p;
                    cache.n = n;
                    pxLightCache.append(cache);
                    if (bounce == settings.GI.traceDepth) /* first hit */
                        visualizeIntensity = cache.i;

                    if (g_debug)
                        printVec3("final gather p:", cache.p);
                }
                bounce--;
            }
            
            rowLightCache.append(pxLightCache);
            renderBuffer.setBlock(x0, y0, step, step, VEC4(ONE, visualizeIntensity));
        }
#pragma omp critical(RenderTask_GI_computeLightCache)
        lightCacheSys_indirect.addLightCache(rowLightCache);
    }
    states.currentRow = rowEnd;

    String progressString = formatProgress(totalRows, rowEnd);
    char buf[512] = { 0 };
    sprintf(buf, "Current [%d]: ", states.resolution_indirectLighting_cur);
    String message = buf;
    message += progressString;
    renderer->updateLog(message);

    if (rowEnd == totalRows)
    {
        states.currentRow = 0;
        states.resolution_indirectLighting_cur++;
        if (states.resolution_indirectLighting_cur > settings.GI.resolution_indirectLighting_max) {
            renderer->logInfo("Finished.");
            sprintf(buf, "Constructing kd tree (%d samples), this may take a while...", lightCacheSys_indirect.numCaches());
            saveLightCachePreview("indirect.png");
            renderer->logInfo(buf);
            timer.tick();
            return true;
        }
        else {
            renderer->logInfo("");
            timer.tick();
            return false;
        }
    }
    else
        return false;
}
bool RenderTask_GI::buildKdTree_indirectLighting(InteractiveRenderer* renderer)
{
    Timer timer;
    timer.tick();
    lightCacheSys_indirect.buildKdTree(settings.GI.knnSamples_indirectLighting, settings.GI.maxKdTreeDepth);
    String dt = formatTime(timer.tick());
    String s = "kd-tree construction finished in ";
    s += dt;
    s += ".";
    renderer->logInfo(s.c_str());
    s = "kd-tree memory consumption: ";
    s += formatSize(lightCacheSys_indirect.computeCacheSize());
    s += ".";
    renderer->logInfo(s.c_str());
    return true;
}
bool RenderTask_GI::visualizeIndirectLighting(InteractiveRenderer* renderer)
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;
    const int& numSamplesPerLight = settings.GI.numSamplesPerLight;

    int row = getNumProcessors();
    int rowStart = states.currentRow;
    int rowEnd = rowStart + row > h ? h : rowStart + row;

    /* uniform scatter in screen space and bounce */
#pragma omp parallel for schedule(static)
    for (int y0 = rowStart; y0 < rowEnd; y0++) {
        XorwowRNG rng(y0, 166);
        Array<light_cache> rowLightCache;
        for (int x0 = 0; x0 < w; x0 ++) {
            if (x0 == 24 && y0 == 582) {
                g_debug = true;
            }
            else
            {
                g_debug = false;
            }
            ray r = camera.raycast(w, h, x0, y0, settings.GI.randomSampling, 1, &rng);
            ray_query Q = scene->rayQuery(&r);
            if (Q.numHits() > 0)
            {
                if (dot(Q[0].normal, r.d) > ZERO)
                    Q[0].normal = -Q[0].normal;

                VEC3 p = Q[0].point + Q[0].normal * REAL(0.001);

                VEC3 avgInten;
                Material* material = scene->getEntityMaterial(Q[0].entityId);
                if ((material->getMaterialType() & Material_Emissive) == false)
                {
                    avgInten = getIndirectLightingIntensity(p, Q[0].normal);
                }
                else {
                    avgInten = material->getEmissiveIntensity(Q[0].uvw, scene).xyz();
                }
                renderBuffer.setPixel(x0, y0, VEC4(ONE, avgInten));
            }
            else {
                renderBuffer.setPixel(x0, y0, settings.render.backgroundColor);
            }
        }
    }
    states.currentRow = rowEnd;

    String progressString = formatProgress(h, rowEnd);
    renderer->updateLog(progressString.c_str());

    if (rowEnd == h) {
        renderer->logInfo("Finished.");
        saveLightCachePreview("indirect_visualization.png");
        states.currentRow = 0;
        /* clear render buffer to prepare for rendering */
        renderBuffer.setBlock(0, 0, renderBuffer.getWidth(), renderBuffer.getHeight(),
            VEC4(ONE, VEC3::all(ZERO)));
        return true;
    }
    else
        return false;
}
bool RenderTask_GI::finalRender(InteractiveRenderer* renderer)
{
    const int& w = settings.render.w;
    const int& h = settings.render.h;
    const int& numSamplesPerLight = settings.GI.numSamplesPerLight;

    int row = getNumProcessors();
    int rowStart = states.currentRow;
    int rowEnd = rowStart + row > h ? h : rowStart + row;

    /* uniform scatter in screen space and bounce */
#pragma omp parallel for schedule(static)
    for (int y0 = rowStart; y0 < rowEnd; y0++) {
        XorwowRNG rng(y0, states.finalRender_curPass);
        Array<light_cache> rowLightCache;
        for (int x0 = 0; x0 < w; x0++) {
            ray r = camera.raycast(w, h, x0, y0, settings.render.antialiasing, 1, &rng);
            ray_query Q;
            int status = findValidRenderPosition(r, Q);
            
            VEC4 finalColor; /* final color for this iteration */

            if (status > 0)
            {
                if (dot(Q[0].normal, r.d) > ZERO)
                    Q[0].normal = -Q[0].normal;
                VEC3 p = Q[0].point + Q[0].normal * REAL(0.001);
                VEC3 n = Q[0].normal;
                Material* material = scene->getEntityMaterial(Q[0].entityId);
                if (material->getMaterialType() & Material_Emissive) {
                    finalColor = material->getEmissiveIntensity(Q[0].uvw, scene);
                }
                else {
                    VEC3 directLighting = scene->computeDirectLightIntensity(p, n, 100, rng);
                    VEC3 indirectLighting = getIndirectLightingIntensity(p, n);
                    VEC4 diffuseColor = material->getDiffuseColor(Q[0].uvw, scene);
                    finalColor = VEC4(ONE, diffuseColor.xyz() * (directLighting + indirectLighting));
                }
            }

            /* average color (Monte Carlo integration) */
            VEC4 prevColor = renderBuffer.getPixel(x0, y0) * REAL(states.finalRender_curPass - 1);
            VEC4 newColor = prevColor + finalColor;
            newColor = newColor / REAL(states.finalRender_curPass);

            renderBuffer.setPixel(x0, y0, newColor);
        }
    }

    states.currentRow = rowEnd;

    char buf[1024];
    sprintf(buf, "Current iteration %d | Elapsed : %s", states.finalRender_curPass,
        formatTime(timer.tick()).c_str());
    renderer->updateLog(buf);

    if (rowEnd == h) {
        states.currentRow = 0;
        states.finalRender_curPass++;
        if (states.finalRender_curPass % 5 == 0) {
            saveLightCachePreview("final_render.png");
        }
    }

    return false;
}
bool RenderTask_GI::checkIfNeedToRefineLighting(int x, int y, int size)
{
    int hsize = (size - 1) / 2;
    
    for (int xx = x - hsize; xx <= x + size; xx++) {
        for (int yy = y - hsize; yy <= y + size; yy++) {
            if (detailBuffer.getPixel(xx, yy).s == ZERO) continue;
            if (detailBuffer.getPixel(xx, yy).x > HALF) return true;
        }
    }
    return false;
}
void RenderTask_GI::computeScreenSpaceDetail(const RenderBuffer& albedoBuffer, 
    const RenderBuffer& normalBuffer, const RenderBuffer& depthBuffer)
{
    const int localSize = 9;/* 4 + 1 + 4 */
    const REAL smoothnessThreshold = REAL(0.98);
    const int& w = settings.render.w;
    const int& h = settings.render.h;

    this->detailBuffer.create(w, h);

    int halfLocalSize = (localSize - 1) / 2;

    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            if (normalBuffer.getPixel(x, y).s == ZERO) continue; /* skip empty pixel */
            
            bool isDetail = false;

            REAL depth = depthBuffer.getPixel(x, y).x;
            VEC3 normal = normalBuffer.getPixel(x, y).xyz();

            /* checking normals in local */
            REAL minDotProduct = ONE;
            for (int xx = x - halfLocalSize; xx <= x + halfLocalSize; xx++) {
                for (int yy = y - halfLocalSize; yy <= y + halfLocalSize; yy++) {
                    if (normalBuffer.getPixel(xx, yy).s == ZERO) {
                        if (xx >= 0 && xx < w && yy >= 0 && yy < h) {
                            isDetail = true;
                        }
                        else continue;
                    }
                    VEC3 n0 = normalBuffer.getPixel(xx, yy).xyz();
                    if (minDotProduct > dot(normal, n0))
                        minDotProduct = dot(normal, n0);
                }
            }
            REAL smoothness = minDotProduct * HALF + HALF;
            if (smoothness < smoothnessThreshold)
                isDetail = true;

            /* checking depth buffer */
            for (int xx = x - 1; xx <= x + 1; xx++) {
                for (int yy = y - 1; yy <= y + 1; yy++) {
                    if (depthBuffer.getPixel(xx, yy).s == ZERO ||
                        normalBuffer.getPixel(xx, yy).s == ZERO)
                    {
                        continue;
                    }
                    REAL d0 = depthBuffer.getPixel(xx, yy).x;
                    VEC3 n0 = normalBuffer.getPixel(xx, yy).xyz();
                    REAL dd = Fabs(d0 - depth);
                    ray r = camera.raycast(w, h, xx, yy);
                    if (dot(r.d, n0) > ZERO)
                        n0 = -n0;
                    REAL theta = Arccos(dot(-r.d, n0));
                    REAL g = camera.convertPixelSize(2, w, depth);
                    REAL dm = g * Tan(theta);
                    if (dd > dm)
                        isDetail = true;
                }
            }

            if (isDetail)
                detailBuffer.setPixel(x, y, VEC4::ones());
            else
                detailBuffer.setPixel(x, y, VEC4(ONE, VEC3::zeros()));
        }
    }
}
VEC3 RenderTask_GI::getIndirectLightingIntensity(const VEC3& p, const VEC3& n)
{
    Array<light_cache> cache = lightCacheSys_indirect.knnSearch(
        p, settings.GI.knnSamples_indirectLighting);
    VEC3 avgInten = VEC3::zeros();
    int acceptedSamples = 0;
    if (g_debug)
        printVec3("p :", p);
    for (int i = 0; i < cache.size(); i++) {
        if (scene->checkVisibilityBetween(p, cache[i].p, n, cache[i].n) == true
            && dot(n, cache[i].n) > REAL(0.7) /* <- only collect coplanar samples*/)
        {
            if (g_debug) 
                printVec3("cache*:", cache[i].i);
            avgInten += cache[i].i;
            acceptedSamples++;
        }
        else
            if (g_debug)
                printVec3("cache :", cache[i].i);
    }
    if (acceptedSamples > 0)
        avgInten /= REAL(acceptedSamples);

    if (g_debug)
        printVec3("inten :", avgInten);

    return avgInten;
}
VEC3 RenderTask_GI::getDirectLightingIntensity(const VEC3& p, const VEC3& n)
{
    Array<light_cache> cache = lightCacheSys_direct.knnSearch(
        p, settings.GI.knnSamples_directLighting);
    VEC3 avgInten = VEC3::zeros();
    int acceptedSamples = 0;
    for (int i = 0; i < cache.size(); i++) {
        if (scene->checkVisibilityBetween(p, cache[i].p, n, cache[i].n) == true
            && dot(n, cache[i].n) > REAL(0.7) /* <- only collect coplanar samples*/)
        {
            avgInten += cache[i].i;
            acceptedSamples++;
        }
    }
    if (acceptedSamples > 0)
        avgInten /= REAL(acceptedSamples);
    return avgInten;
}
int RenderTask_GI::findValidLightCachePosition(ray& r, ray_query& Q)
{
    int maxTries = 16;
    while (maxTries > 0) { 
        /* here i can write an infinite loop but for 
        robustness i still setup a max iteration value */
        
        Q = scene->rayQueryIgnoreMaterial(&r, Material_Volumetric);
        if (Q.numHits() > 0) {
            Material* material = scene->getEntityMaterial(Q[0].entityId);
            unsigned int materialType = material->getMaterialType();
            if ((materialType & Material_Emissive) || !(materialType & Material_Diffuse)) {
                return 0; /* hit a surface but it's emissive or not diffusive */
            }
            else {
                /* a diffusive surface, check its opacity */
                VEC4 albedo = material->getDiffuseColor(Q[0].uvw, scene);
                if (albedo.s < HALF) {
                    /* leak through */
                    if (dot(r.d, Q[0].normal) < ZERO) {
                        r.o = Q[0].point - Q[0].normal * REAL(0.001);
                    }
                    else {
                        r.o = Q[0].point + Q[0].normal * REAL(0.001);
                    }
                }
                else {
                    return 1; /* we find a valid cache point */
                }
            }
        }
        else
            return -1;
        maxTries--;
    }
    return -1;
}
int RenderTask_GI::findValidRenderPosition(ray& r, ray_query& Q)
{
    int maxTries = 64;
    while (maxTries > 0) {
        /* here i can write an infinite loop but for
        robustness i still setup a max iteration value */
        Q = scene->rayQuery(&r);
        if (Q.numHits() > 0) {
            Material* material = scene->getEntityMaterial(Q[0].entityId);
            unsigned int materialType = material->getMaterialType();
            
            if (materialType & Material_Diffuse) {
                /* a diffusive surface, check its opacity */
                VEC4 albedo = material->getDiffuseColor(Q[0].uvw, scene);
                if (albedo.s < HALF) {
                    /* leak through */
                    if (dot(r.d, Q[0].normal) < ZERO) {
                        r.o = Q[0].point - Q[0].normal * REAL(0.001);
                    }
                    else {
                        r.o = Q[0].point + Q[0].normal * REAL(0.001);
                    }
                }
                else {
                    return 1; /* we find a valid cache point */
                }
            }
            else {
                /* hit other types of surface */
                return 1;
            }
        }
        else
            return -1;
        maxTries--;
    }
    return -1;
}
void RenderTask_GI::createRenderTask(
    Scene * scene, Camera * camera, RenderSettings* settings,
    RenderTask_Albedo* albedoTask,
    RenderTask_Normal* normalTask,
    RenderTask_Depth* depthTask)
{
    this->displayTexture.gen(settings->render.w, settings->render.h, GL_RGBA, GL_NEAREST);
    this->scene = scene;
    this->albedoTask = albedoTask;
    this->normalTask = normalTask;
    this->depthTask = depthTask;
    this->settings = *settings;
    this->camera = *camera;
    this->renderBuffer.create(settings->render.w, settings->render.h);
    this->mainRng.init(1732, 2906);
    this->settings.validateParameters();
    states.currentRow = 0;
    states.isDirectLightingComputed = false;
    states.isDirectLightingKdTreeBuilt = false;
    states.isIndirectLightingComputed = false;
    states.isIndirectLightingKdTreeBuilt = false;
    states.isIndirectLightingVisualized = false;
    states.isFinalRenderCompleted = false;
    states.resolution_indirectLighting_cur = this->settings.GI.resolution_indirectLighting_min;
    states.resolution_directLighting_cur = this->settings.GI.resolution_directLighting_min;
    states.isScreenSpaceDetailComputed = false;
    states.finalRender_curPass = 1;
}
void RenderTask_GI::interactiveSubmit(InteractiveRenderer * renderer)
{
    /* submit rendered result to GPU texture, so that it can be seen */
    /* from the users. */
    const BYTE* p = this->renderBuffer.asBitmap(colorToRGBA);
    int stride = sizeof(BYTE) * 4 * this->renderBuffer.getWidth();
    int w = this->renderBuffer.getWidth();
    int h = this->renderBuffer.getHeight();
    int bytes = stride * h;
    BYTE* pc = (BYTE*)malloc(bytes);
    memcpy(pc, p, bytes);
    /* draw annotations */
    if (states.isDirectLightingComputed == false) {
        const int step = ipow(2,-states.resolution_directLighting_cur);
        memset(pc + states.currentRow * step  * stride, 255, stride);
    }
    else if (states.isIndirectLightingComputed == false) {
        const int step = ipow(2, -states.resolution_indirectLighting_cur);
        memset(pc + states.currentRow * step * stride, 255, stride);
    }
    else if (states.isIndirectLightingVisualized == false) {
        memset(pc + states.currentRow * stride, 255, stride);
    }
    /*else if (states.isFinalRenderCompleted == false) {
        memset(pc + states.currentRow * stride, 255, stride);
    }*/
    /* GL texture is vertically flipped, so we need to flip them too */
    vflipBuffer(pc, stride, h);
    displayTexture.update(pc);
    free(pc);
}
RenderStatus RenderTask_GI::interactiveRender(InteractiveRenderer * renderer)
{
    if (states.isScreenSpaceDetailComputed == false) {
        this->computeScreenSpaceDetail(
            albedoTask->getRenderBuffer(),
            normalTask->getRenderBuffer(),
            depthTask->getRenderBuffer());
        //detailBuffer.saveBitmap(colorToRGBA, "detail.png");
        states.isScreenSpaceDetailComputed = true;
        return Render_Running;
    }
    else if (states.isDirectLightingComputed == false) {
        states.isDirectLightingComputed = computeDirectLighting(renderer);
        return Render_Running;
    }
    else if (states.isDirectLightingKdTreeBuilt == false) {
        states.isDirectLightingKdTreeBuilt = buildKdTree_directLighting(renderer);
        return Render_Running;
    }
    else if (states.isIndirectLightingComputed == false) {
        states.isIndirectLightingComputed = computeIndirectLighting(renderer);
        return Render_Running;
    }
    else if (states.isIndirectLightingKdTreeBuilt == false) {
        states.isIndirectLightingKdTreeBuilt = buildKdTree_indirectLighting(renderer);
        return Render_Running;
    }
    else if (states.isIndirectLightingVisualized == false) {
        states.isIndirectLightingVisualized = visualizeIndirectLighting(renderer);
        return Render_Running;
    }
    else if (states.isFinalRenderCompleted == false) {
        states.isFinalRenderCompleted = finalRender(renderer);
        return Render_Running;
    }
    /*
    else if (xxx=false){}
    */
    return Render_Finished;
}
void RenderTask_GI::onRenderStart(InteractiveRenderer * renderer)
{
    RenderTask::onRenderStart(renderer);
    renderer->logInfo("Start GI rendering...");
    renderer->logInfo("Computing light cache...");
    renderer->logInfo("");
}
void RenderTask_GI::onRenderFinish(InteractiveRenderer * renderer)
{
    RenderTask::onRenderFinish(renderer);

    renderer->logInfo("Render finished.");
    double dt = timer.tick();
    char buffer[512];
    sprintf(buffer, "Render time: %.2lfsec.", dt);
    renderer->logInfo(buffer);
}
bool RenderTask_GI::checkRenderTask(String & log)
{
    if (scene->getEmissiveEntityIds().size() == 0) {
        log = "No light in scene.";
        return false;
    }
    return true;
}

void LightCacheSystem::addLightCache(const Array<light_cache>& cache)
{
    lightCache.append(cache);
}
void LightCacheSystem::buildKdTree(const int & k, const int & maxDepth)
{
    destroyKdTree();
    this->kdSplitThreshold = k * 2;
    this->kdMaxDepth = maxDepth;

    kdroot = new LightCacheKdTreeNode();
    kdroot->depth = 1;
    buildKdTree_boxFit(0, lightCache.size(), &(kdroot->box.bmin), &(kdroot->box.bmax));
    buildKdTree_recur(kdroot, 0, lightCache.size());
    printf(" complete.\n");
}
bool LightCacheSystem::save(const char * file) const
{
    FILE* fp = NULL;
    if ((fp = fopen(file, "wb")) == NULL)
        return false;

    const char magicString[64] = "rkJYme3skGz11EQtpwTuCh";
    const char className[64] = "LightCacheSystem";
    fwrite(magicString, 1, 64, fp);
    fwrite(className, 1, 64, fp);

    int iBuffer;
    float fBuffer;

    /* write all vertices (save as single precision float) */
    iBuffer = lightCache.size(); fwrite(&iBuffer, 4, 1, fp);
    for (int i = 0; i < lightCache.size(); i++) {
        VEC3 data;
        data = lightCache[i].p;
        fBuffer = float(data.x); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.y); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.z); fwrite(&fBuffer, 4, 1, fp);
        data = lightCache[i].n;
        fBuffer = float(data.x); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.y); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.z); fwrite(&fBuffer, 4, 1, fp);
        data = lightCache[i].i;
        fBuffer = float(data.x); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.y); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.z); fwrite(&fBuffer, 4, 1, fp);
    }

    fclose(fp);
    return true;
}
bool LightCacheSystem::load(const char * file)
{
    lightCache.clear();

    FILE* fp = NULL;
    if ((fp = fopen(file, "rb")) == NULL)
        return false;

    const char magicString[64] = "rkJYme3skGz11EQtpwTuCh";
    const char className[64] = "LightCacheSystem";

    char buf[64];
    fread(buf, 1, 64, fp);
    if (strcmp(buf, magicString) != 0) {
        fclose(fp);
        return false;
    }
    fread(buf, 1, 64, fp);
    if (strcmp(buf, className) != 0) {
        fclose(fp);
        return false;
    }

    int cacheSize;

    /* read all (save as single precision float) */
    fread(&cacheSize, 4, 1, fp);
    lightCache.reserve(cacheSize);
    for (int i = 0; i < cacheSize; i++) {
        light_cache cache;
        fread(&cache.p.x, 4, 1, fp);
        fread(&cache.p.y, 4, 1, fp);
        fread(&cache.p.z, 4, 1, fp);
        fread(&cache.n.x, 4, 1, fp);
        fread(&cache.n.y, 4, 1, fp);
        fread(&cache.n.z, 4, 1, fp);
        fread(&cache.i.x, 4, 1, fp);
        fread(&cache.i.y, 4, 1, fp);
        fread(&cache.i.z, 4, 1, fp);
        lightCache.append(cache);
    }

    fclose(fp);
    return true;
}

Array<light_cache> LightCacheSystem::knnSearch(const VEC3& query, 
    const int& k, bool approximate) const
{
    /* find initial solution for knn */
    PriorityQueue<int, REAL> solution;
    LightCacheKdTreeNode* searchedLeaf = NULL;
    knn_findInitialSolution(query, k, solution, &searchedLeaf);

    /* if performing an approximate k-NN search, just use initial solution */
    /* as our final solution */

    if (approximate == false) {
        /* else we refine the solution */
        /* build hypersphere for refining knn solution */
        int iTop; REAL dTop;
        solution.query(iTop, dTop);
        sphere hyperSphere;
        hyperSphere.r = Sqrt(dTop);
        hyperSphere.p = query;

        knn_refineSolution(hyperSphere, k, searchedLeaf, solution);
    }

    Array<light_cache> cache;
    /* collect light cache points and return results */
    for (int i = 0; i < solution.size(); i++) {
        /* note that here we use solution.size() instead of k because */
        /* the actual solution can contain number of cache points smaller than k */
        int iCache = solution[i];
        cache.append(lightCache[iCache]);
    }

    return cache;
}
int LightCacheSystem::numCaches() const
{
    return lightCache.size();
}
int LightCacheSystem::computeCacheSize() const
{
    /* return memory consumption in bytes */
    int totalMem = 0;

    if (kdroot == NULL) return 0;

    Stack<LightCacheKdTreeNode*> stack;
    stack.push(kdroot);
    while (stack.isEmpty() == false)
    {
        LightCacheKdTreeNode* topNode;
        stack.pop(topNode);
        totalMem += sizeof(LightCacheKdTreeNode);
        totalMem += sizeof(int) * topNode->n;
        if (topNode->left)
            stack.push(topNode->left);
        if (topNode->right)
            stack.push(topNode->right);
    }

    return totalMem;
}
int LightCacheSystem::maxNodeCacheCount() const
{
    /* return memory consumption in bytes */
    int maxCount = 0;

    if (kdroot == NULL) return 0;

    Stack<LightCacheKdTreeNode*> stack;
    stack.push(kdroot);
    while (stack.isEmpty() == false)
    {
        LightCacheKdTreeNode* topNode;
        stack.pop(topNode);
        if (maxCount < topNode->n)
            maxCount = topNode->n;
        if (topNode->left)
            stack.push(topNode->left);
        if (topNode->right)
            stack.push(topNode->right);
    }

    return maxCount;
}
int LightCacheSystem::minNodeCacheCount() const
{
    /* return memory consumption in bytes */
    int minCount = numCaches();

    if (kdroot == NULL) return 0;

    Stack<LightCacheKdTreeNode*> stack;
    stack.push(kdroot);
    while (stack.isEmpty() == false)
    {
        LightCacheKdTreeNode* topNode;
        stack.pop(topNode);
        if (topNode->left == NULL && topNode->right == NULL && minCount > topNode->n)
            minCount = topNode->n;
        if (topNode->left)
            stack.push(topNode->left);
        if (topNode->right)
            stack.push(topNode->right);
    }

    return minCount;
}

int _lightCacheQuickSort_pivot(Array<light_cache>& arr, const int& L, const int& R, 
    const int& axis)
{
    /* sort an array from small to large */
    if (R - L <= 1) return L;
    const int pivot = R - 1;
    int l = R - 1, r = -1;
    for (int i = L; i < R - 1; i++) {
        if (arr[i].p.e[axis] > arr[pivot].p.e[axis]) {
            /* if cur > pivot, set left */
            l = i;
            /* find an element smaller than pivot */
            bool find = false;
            for (int j = i + 1; j < R - 1; j++) {
                r = j;
                if (arr[r].p.e[axis] < arr[pivot].p.e[axis]) {
                    swap(arr[l], arr[r]);
                    find = true;
                    break;
                }
            }
            if (!find) {
                /* all elements between left and pivot are larger than pivot */
                /* just break the loop to swap left and pivot and we are done. */
                break;
            }
        }
    }
    swap(arr[l], arr[pivot]);
    /* now all elements smaller than the pivot are on the left side, */
    /* all elements larger than the pivot are on the right side */
    return l;
}
void _lightCacheQuickSort(Array<light_cache>& arr, int start, int end, int axis)
{
    if (start < end) {
        int* stack = new int[2 * (end - start + 1)];
        int top = -1;

        stack[++top] = start;
        stack[++top] = end;

        while (top >= 0) {
            end = stack[top--];
            start = stack[top--];

            int mid = _lightCacheQuickSort_pivot(arr, start, end, axis);
            if (mid - start > 1) {
                /* still contain some elements on the left */
                stack[++top] = start;
                stack[++top] = mid;
            }
            if (end - mid > 1) {
                /* still contain some elements on the right */
                stack[++top] = mid + 1;
                stack[++top] = end;
            }
        }
        delete[] stack;
    }
}

void LightCacheSystem::sortLightCache(int start, int count, int axis)
{
    _lightCacheQuickSort(this->lightCache, start, start + count, axis);
}
void LightCacheSystem::destroyKdTree_recur(LightCacheKdTreeNode * node)
{
    if (node == NULL)
        return;
    if (node->left)
    {
        destroyKdTree_recur(node->left);
        node->left = NULL;
    }
    if (node->right) {
        destroyKdTree_recur(node->right);
        node->right = NULL;
    }
    if (node->idxs) {
        free(node->idxs);
        node->idxs = NULL;
    }
    delete node;
}
void LightCacheSystem::destroyKdTree()
{
    destroyKdTree_recur(kdroot);
    kdroot = NULL;
    _buildCounter = 0;
}
void LightCacheSystem::buildKdTree_splitBox(const aabb * box, const char axis, const REAL & value, aabb * left, aabb * right)
{
    switch (axis) {
    case 0:
        left->bmin = box->bmin;
        left->bmax = VEC3(value, box->bmax.y, box->bmax.z);
        right->bmin = VEC3(value, box->bmin.y, box->bmin.z);
        right->bmax = box->bmax;
        break;
    case 1:
        left->bmin = box->bmin;
        left->bmax = VEC3(box->bmax.x, value, box->bmax.z);
        right->bmin = VEC3(box->bmin.x, value, box->bmin.z);
        right->bmax = box->bmax;
        break;
    case 2:
        left->bmin = box->bmin;
        left->bmax = VEC3(box->bmax.x, box->bmax.y, value);
        right->bmin = VEC3(box->bmin.x, box->bmin.y, value);
        right->bmax = box->bmax;
        break;
    default:
        break;
    }
}
void LightCacheSystem::buildKdTree_boxFit(const int& start, const int& count, VEC3 * bmin, VEC3 * bmax)
{
    VEC3 _bmin = VEC3(REAL(REAL_MAX), REAL(REAL_MAX), REAL(REAL_MAX));
    VEC3 _bmax = VEC3(REAL(-REAL_MAX), REAL(-REAL_MAX), REAL(-REAL_MAX));

    for (int i = start; i < start+count; i++) {
        VEC3 p = lightCache[i].p;
        if (_bmin.x > p.x) _bmin.x = p.x;
        if (_bmin.y > p.y) _bmin.y = p.y;
        if (_bmin.z > p.z) _bmin.z = p.z;
        if (_bmax.x < p.x) _bmax.x = p.x;
        if (_bmax.y < p.y) _bmax.y = p.y;
        if (_bmax.z < p.z) _bmax.z = p.z;
    }

    *bmin = _bmin;
    *bmax = _bmax;
}
void LightCacheSystem::buildKdTree_recur(LightCacheKdTreeNode* node, 
    const int& start, const int& count)
{
    bool shouldSplit = true;
    aabb lBox, rBox;
    int leftCount, leftStart, leftEnd;
    int rightCount, rightStart, rightEnd;

    const REAL delta = REAL(0.01);

    /* find dimension with largest span */
    VEC3 dxyz = node->box.bmax - node->box.bmin;
    int axis = argmax3(dxyz.x, dxyz.y, dxyz.z);

    /* sort array within axis with largest span */
    this->sortLightCache(start, count, axis);
    
    /* determine if this node should split */
    if (node->depth < this->kdMaxDepth && count > this->kdSplitThreshold && dxyz.e[axis] > delta)
    {
        shouldSplit = true;

        /* since the array is already sorted, we just need to find previous n/2 points */
        leftCount = (count + 1) / 2;
        rightCount = count - leftCount;
        leftStart = start;
        leftEnd = start + leftCount;
        rightStart = leftEnd;
        rightEnd = rightStart + rightCount;

        REAL midvalue = HALF * (lightCache[rightStart].p.e[axis] + lightCache[rightStart - 1].p.e[axis]);

        /* split box and put samples into different sub-boxes */
        buildKdTree_splitBox(&(node->box), axis, midvalue, &lBox, &rBox);
        
    }
    else {
        shouldSplit = false;
    }

    /* * * split topNode into two child nodes or turn it to leaf topNode * * */

    if (shouldSplit) {
        /* build nodes */
        node->left = new LightCacheKdTreeNode();
        node->right = new LightCacheKdTreeNode();
        node->left->depth = node->depth + 1;
        node->right->depth = node->depth + 1;
        node->left->box = lBox;
        node->right->box = rBox;

        /* resursive build */
        buildKdTree_recur(node->left, leftStart, leftCount);
        buildKdTree_recur(node->right, rightStart, rightCount);
    }
    else {
        node->n = count;
        node->idxs = (int*)malloc(sizeof(int)*(node->n));
        for (int i = 0; i < count; i++) {
            node->idxs[i] = start + i;
        }
        _buildCounter += count;
    }
    printf("\rbuilding cache [%d/%d]...", _buildCounter, this->numCaches());
}
void LightCacheSystem::knn_findInitialSolution(const VEC3 & query, 
    const int& k, PriorityQueue<int, REAL>& solution, LightCacheKdTreeNode** searchedLeaf) const
{
    *searchedLeaf = NULL;

    /* find nearest leaf from query point */
    VEC3 pBox;
    REAL dist = squaredDistanceBetween(query, kdroot->box, pBox);
    if (dist == ZERO)
        /* query point is inside the bounding box */
        pBox = query;

    /* now find all leaf nodes that pBox lies in, */
    /* and among these leaves, find nearest k points */
    /* as our initial solution for k-NN search */

    Stack<LightCacheKdTreeNode*> nodeStack;
    nodeStack.reserve(kdMaxDepth + 1);
    nodeStack.push(kdroot); /* initialize stack and trace */

    while (!nodeStack.isEmpty()) {
        LightCacheKdTreeNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                int iCache = topNode->idxs[i];
                VEC3 pCache = lightCache[iCache].p;
                REAL dCache = lensq(pBox - pCache);
                if (solution.size() < k)
                    /* if priority queue not full */
                    /* just add this point as part */
                    /* of the solution */
                    solution.put(iCache, dCache);
                else{
                    /* check if current cache point */
                    /* is nearer than current furthest cache point */
                    /* if so, update solution */
                    int iFurthest;
                    REAL dFurthest;
                    solution.query(iFurthest, dFurthest);
                    if (dCache < dFurthest)
                        solution.replace(iCache, dCache);
                }
            }
            /* note that for solving the initial solution,  */
            /* we only need to check one leaf node, now the */
            /* leaf has been checked, we break the outer    */
            /* while loop to end the search process.        */
            *searchedLeaf = topNode;
            break;
        }
        else {
            /* non leaf node */
            if (topNode->left && isAinB(pBox, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isAinB(pBox, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
}
void LightCacheSystem::knn_refineSolution(const sphere & sph, const int& k, 
    const LightCacheKdTreeNode* searchedLeaf, PriorityQueue<int, REAL>& solution) const
{
    Stack<LightCacheKdTreeNode*> nodeStack;
    nodeStack.reserve(kdMaxDepth + 1);
    nodeStack.push(kdroot); /* initialize stack and trace */

    while (!nodeStack.isEmpty()) {
        LightCacheKdTreeNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            if (topNode == searchedLeaf) {
                /* if this leaf node is already searched in */
                /* initial solution stage, ignore it */
                continue;
            }
            for (int i = 0; i < topNode->n; i++) {
                int iCache = topNode->idxs[i];
                VEC3 pCache = lightCache[iCache].p;
                REAL dCache = lensq(sph.p - pCache);
                /* check if current cache point */
                /* is nearer than current furthest cache point */
                /* if so, update solution */
                int iFurthest;
                REAL dFurthest;
                solution.query(iFurthest, dFurthest);
                if (dCache < dFurthest)
                    solution.replace(iCache, dCache);
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(sph, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(sph, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
}
const light_cache& LightCacheSystem::getCache(const int& index) const
{
    return lightCache[index];
}
LightCacheSystem::LightCacheSystem()
{
    kdroot = NULL;
    this->kdMaxDepth = 0;
    this->kdSplitThreshold = 0;
    this->kMax = 0;
    _buildCounter = 0;
}
LightCacheSystem::~LightCacheSystem()
{
    destroyKdTree();
}
void LightCacheSystem::addLightCache(const light_cache & cache)
{
    lightCache.append(cache);
}
LightCacheSystem::LightCacheKdTreeNode::LightCacheKdTreeNode()
{
    left = right = NULL;
    depth = n = 0;
    idxs = NULL;
}

ray Camera::raycast(const int& w, const int& h, const int& x, const int& y, 
    bool randomSampling, const int& subpixels, XorwowRNG* rng)
{
    /* calculate raytrace parameters */
    VEC3 F = normalize(look - pos);   /* front */
    VEC3 R = normalize(cross(F, up)); /* right */
    VEC3 U = normalize(cross(R, F));  /* up */
    VEC3 X = 2 * Tan(fov / TWO) * R;
    VEC3 Y = len(X) * REAL(h) / REAL(w) * (-U);
    VEC3 O = pos + F - X / TWO - Y / TWO;
    VEC3 dX = X / REAL(w);
    VEC3 dY = Y / REAL(h);

    REAL j = REAL(subpixels) * HALF;
    REAL jx = ZERO, jy = ZERO;
    if (randomSampling) {
        jx = rng->uniform(-j, j);
        jy = rng->uniform(-j, j);
    }
    VEC3 T = O + dX * (REAL(x) + j + jx) + dY * (REAL(y) + j + jy);
    ray r(pos, T - pos);

    return r;
}

REAL Camera::convertPixelSize(int px, int w, REAL viewDistance)
{
    REAL radPerPx = fov / REAL(w);
    return radPerPx * REAL(px) * viewDistance;
}

void Camera::worldPositionToScreen(VEC3 p, int w, int h, int& x, int& y)
{
    /* calculate raytrace parameters */
    VEC3 F = normalize(look - pos);   /* front */
    VEC3 R = normalize(cross(F, up)); /* right */
    VEC3 U = normalize(cross(R, F));  /* up */
    VEC3 X = 2 * Tan(fov / TWO) * R;
    VEC3 Y = len(X) * REAL(h) / REAL(w) * (-U);
    VEC3 O = pos + F - X / TWO - Y / TWO;
    VEC3 dX = X / REAL(w);
    VEC3 dY = Y / REAL(h);

    VEC3 S = pos - p;
    if (len(S) < REAL(0.001)) {
        /* too near */
        x = y = -1;
        return;
    }
    REAL dotThresh = dot(F, normalize(O - pos));
    if (dot(F, normalize(-S)) < dotThresh) {
        /* outside screen */
        x = y = -1;
        return;
    }

    /* scale S to view plane */
    REAL t = len(VecProject(S, F));
    S /= t;

    /* project S to R */
    VEC3 A = VecProject(S, R);
    VEC3 B = VecProject(S, -U);
    VEC3 dA = X / TWO - A;
    VEC3 dB = Y / TWO - B;
    if (dot(dA, X / TWO) < ZERO || dot(dB, Y / TWO) < ZERO)
    {
        x = y = -1;
        return;
    }
    x = int(len(dA) / len(dX));
    y = int(len(dB) / len(dY));
}
