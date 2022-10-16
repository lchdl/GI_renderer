#pragma once
#include "basedefs.h"
#include "basemath.h"
#include "geometry/tribox.h" /* triangle AABB intersect test */
#include "geometry/tritri.h" /* triangle intersect test */
#include "random.h"

typedef VEC3 point;
struct sphere { VEC3 p; REAL r; };
struct segment { VEC3 p[2]; };
/* axis aligned bounding box */
struct aabb { VEC3 bmin, bmax; };             
struct triangle { VEC3 p[3]; };
/* similar with triangle but with more information */
struct triangleEx { 
    VEC3 p[3], n[3], t[3]; 
    operator triangle() const;
}; 
struct plane { VEC3 p, n; };
struct ray { 
    VEC3 o, d; 
    ray();
    ray(VEC3 o, VEC3 d);
};
struct ray_hit {
    bool hit;
    int entityId;
    int triangleId;
    VEC3 point, uvw, normal;
    REAL dist;
};

/* generate random points on a sphere/hemisphere */
VEC3 randHemi();
VEC3 randHemi(XorwowRNG & rng);
VEC3 randHemi(VEC3 n);
VEC3 randHemi(VEC3 n, XorwowRNG & rng);
VEC3 randSphere();
VEC3 randSphere(XorwowRNG & rng);

bool isIntersect(const triangle& A, const aabb& B);                           /* <- fast & robust! */
bool isIntersect(const aabb& A, const aabb& B);
bool isIntersect(const ray& A, const aabb& B);                                /* <- fast & robust! */
bool isIntersect(const ray& A, const triangle& B);                            /* <- fast & robust! */
bool isIntersect(const ray& A, const triangle& B, REAL& u, REAL& v, REAL& t); /* <- fast & robust! */
bool isIntersect(const triangle& A, const triangle& B);
bool isIntersect(const triangle& A, const triangle& B, int& isCoplanar, segment& s);
bool isIntersect(const sphere& A, const triangle& B);
bool isIntersect(const triangle& A, const plane& B);
bool isIntersect(const sphere& A, const aabb& B);

aabb computeBoundingBoxFromPointCloud(const Array<VEC3>& points);

REAL distanceBetween(const sphere& A, const aabb& B);
REAL distanceBetween(const point& A, const aabb& B);
REAL distanceBetween(const point& A, const triangle& B);
REAL distanceBetween(const point& A, const aabb& B, VEC3& p);
REAL squaredDistanceBetween(const point & A, const triangle & B);
REAL squaredDistanceBetween(const point & A, const aabb & B, VEC3 & p);
REAL squaredDistanceBetween(const point & A, const aabb & B);
REAL surfaceArea(const aabb& A);
REAL volume(const aabb& A);
bool isAinB(const point& A, const aabb& B);
bool segIntersect1d(const REAL& x1, const REAL& x2, const REAL& y1, const REAL& y2, REAL& r1, REAL& r2);
bool aabbIntersect(const aabb& A, const aabb& B, aabb& result);
/*
-1: all vertices smaller than value
0: triangle intersects with plane
+1: all vertices larger than value
*/
int triangleAxisPlaneTest(const triangle& A, const char& axis, const REAL& value); 


