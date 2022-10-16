#include "geometry.h"

/* generate random points on a hemisphere */
VEC3 randHemi()
{
    VEC3 d;
    REAL l;
    do {
        d.x = uniform(REAL(-1.0), REAL(+1.0));
        d.y = uniform(REAL(-1.0), REAL(+1.0));
        d.z = uniform(REAL(-1.0), REAL(+1.0));
        l = lensq(d);
    } while (l > ONE || l < REAL(0.0001)); /* reject sampling */
    d.z = Fabs(d.z);
    return normalize(d);
}

VEC3 randSphere()
{
    VEC3 d;
    REAL l;
    do {
        d.x = uniform(REAL(-1.0), REAL(+1.0));
        d.y = uniform(REAL(-1.0), REAL(+1.0));
        d.z = uniform(REAL(-1.0), REAL(+1.0));
        l = lensq(d);
    } while (l > ONE || l < REAL(0.0001)); /* reject sampling */
    return normalize(d);
}

VEC3 randSphere(XorwowRNG & rng)
{
    VEC3 d;
    REAL l;
    do {
        d.x = rng.uniform(REAL(-1.0), REAL(+1.0));
        d.y = rng.uniform(REAL(-1.0), REAL(+1.0));
        d.z = rng.uniform(REAL(-1.0), REAL(+1.0));
        l = lensq(d);
    } while (l > ONE || l < REAL(0.0001)); /* reject sampling */
    return normalize(d);
}

VEC3 randHemi(XorwowRNG& rng)
{
    VEC3 d;
    REAL l;
    do {
        d.x = rng.uniform(REAL(-1.0), REAL(+1.0));
        d.y = rng.uniform(REAL(-1.0), REAL(+1.0));
        d.z = rng.uniform(REAL(-1.0), REAL(+1.0));
        l = lensq(d);
    } while (l > ONE || l < REAL(0.0001)); /* reject sampling */
    d.z = Fabs(d.z);
    return normalize(d);
}

VEC3 randHemi(VEC3 n)
{
    /* build local coordinate axis */
    VEC3 u = VEC3(ZERO, ZERO, ONE);
    VEC3 t = cross(n, u);
    if (lensq(t) < REAL(0.001)) { u = VEC3(ZERO, ONE, ZERO); }
    VEC3 v = normalize(cross(n, u));
    u = normalize(cross(v, n));
    /* sample from hemisphere */
    t = randHemi();
    // transform t to new coordinate
    VEC3 d;
    d.x = u.x * t.x + v.x * t.y + n.x * t.z;
    d.y = u.y * t.x + v.y * t.y + n.y * t.z;
    d.z = u.z * t.x + v.z * t.y + n.z * t.z;
    return d;
}

VEC3 randHemi(VEC3 n, XorwowRNG& rng)
{
    /* build local coordinate axis */
    VEC3 u = VEC3(ZERO, ZERO, ONE);
    VEC3 t = cross(n, u);
    if (lensq(t) < REAL(0.001)) { u = VEC3(ZERO, ONE, ZERO); }
    VEC3 v = normalize(cross(n, u));
    u = normalize(cross(v, n));
    /* sample from hemisphere */
    t = randHemi(rng);
    // transform t to new coordinate
    VEC3 d;
    d.x = u.x * t.x + v.x * t.y + n.x * t.z;
    d.y = u.y * t.x + v.y * t.y + n.y * t.z;
    d.z = u.z * t.x + v.z * t.y + n.z * t.z;
    return d;
}

REAL distanceBetween(const sphere & A, const aabb & B)
{
    return distanceBetween(A.p, B) - A.r;
}

REAL distanceBetween(const point & A, const aabb & B)
{
    return Sqrt(squaredDistanceBetween(A, B));
}

REAL distanceBetween(const point & A, const triangle & B)
{
    return Sqrt(squaredDistanceBetween(A, B));
}

REAL distanceBetween(const point & A, const aabb & B, VEC3 & p)
{
    return Sqrt(squaredDistanceBetween(A, B, p));
}

void _calcTriangleNormals(
    /* in */
    const VEC3& v1, const VEC3& v2, const VEC3& v3,
    /* out */
    VEC3& n, VEC3& n1, VEC3& n2, VEC3& n3
) {
    VEC3 v1v2 = v2 - v1;
    VEC3 v1v3 = v3 - v1;
    VEC3 v2v3 = v3 - v2;
    n = normalize(cross(v1v2, v1v3));
    n1 = normalize(cross(v1v2, n));
    n2 = normalize(cross(v2v3, n));
    n3 = normalize(cross(n, v1v3));
}

void _projectPointToPlane(
    /* in */
    const VEC3 & p, const VEC3 & o, const VEC3 & n,
    /* out */
    VEC3& j, int& s) {

    VEC3 po = o - p;
    REAL projlen = dot(n, po);
    VEC3 pj = n * projlen;
    j = p + pj;
    s = (projlen < 0) ? (+1) : (-1);
}

REAL squaredDistanceBetween(const point & A, const triangle & B)
{
    const VEC3& v1 = B.p[0];
    const VEC3& v2 = B.p[1];
    const VEC3& v3 = B.p[2];

    VEC3 n, n1, n2, n3;
    _calcTriangleNormals(v1, v2, v3, n, n1, n2, n3);

    /* project point to triangle */
    VEC3 pj;
    int _s = 0; /* unused */
    _projectPointToPlane(A, v1, n, pj, _s);

    /* project point to edge planes */
    VEC3 j[3];
    int s[3];
    _projectPointToPlane(pj, v1, n1, j[0], s[0]);
    _projectPointToPlane(pj, v2, n2, j[1], s[1]);
    _projectPointToPlane(pj, v3, n3, j[2], s[2]);

    VEC3 q;

    if (s[0] < 0 && s[1] < 0 && s[2] < 0) {
        /* projected point j is in the triangle */
        q = pj;
    }
    else {
        int S = s[0] * s[1] * s[2];
        if (S == -1) {
            /* nearest point is one of the three vertices of the triangle */
            REAL d0 = REAL_MAX;
            int i0;
            for (int i = 0; i < 3; i++) {
                REAL d1 = lensq(pj - B.p[i]);
                if (d0 > d1) {
                    d0 = d1;
                    i0 = i;
                }
            }
            q = B.p[i0];
        }
        else { /* S == +1 */
            /* nearest point is at triangle vertices or edges */
            for (int i = 0; i < 3; i++) {
                if (s[i] < 0) continue; /* this will ignore two edges */
                VEC3 vj = j[i] - B.p[i];
                VEC3 vv = B.p[(i + 1) % 3] - B.p[i];
                REAL a = dot(vv, vv);
                REAL b = dot(vj, vv);
                if (b < 0) { /* B<0, v[i] is the nearest point */
                    q = B.p[i];
                }
                else if (b < a) { /* 0<B<A, j[i] is the nearest point */
                    q = j[i];
                }
                else { /* B>A, v[(i + 1) % 3] is the nearest point */
                    q = B.p[(i + 1) % 3];
                }
            }
        }
    }
    return lensq(A - q);
}

REAL squaredDistanceBetween(const point & A, const aabb & B, VEC3& p)
{
    VEC3 xChoices = VEC3(B.bmin.x, A.x, B.bmax.x);
    VEC3 yChoices = VEC3(B.bmin.y, A.y, B.bmax.y);
    VEC3 zChoices = VEC3(B.bmin.z, A.z, B.bmax.z);
    p = VEC3(
        xChoices.e[argmax3(B.bmin.x - A.x, ZERO, A.x - B.bmax.x)],
        yChoices.e[argmax3(B.bmin.y - A.y, ZERO, A.y - B.bmax.y)],
        zChoices.e[argmax3(B.bmin.z - A.z, ZERO, A.z - B.bmax.z)]);
    /* p is the nearest point on the box */
    return lensq(A - p);
}

REAL squaredDistanceBetween(const point & A, const aabb & B)
{
    VEC3 p;
    return squaredDistanceBetween(A, B, p);
}

REAL surfaceArea(const aabb & A)
{
    REAL dx = A.bmax.x - A.bmin.x;
    REAL dy = A.bmax.y - A.bmin.y;
    REAL dz = A.bmax.z - A.bmin.z;
    return REAL(2.0) * (dx * dy + dy * dz + dz * dx);
}

REAL volume(const aabb & A)
{
    REAL dx = A.bmax.x - A.bmin.x;
    REAL dy = A.bmax.y - A.bmin.y;
    REAL dz = A.bmax.z - A.bmin.z;
    return dx * dy * dz;
}

bool isAinB(const point & A, const aabb & B)
{
    if (B.bmin.x <= A.x && A.x <= B.bmax.x &&
        B.bmin.y <= A.y && A.y <= B.bmax.y &&
        B.bmin.z <= A.z && A.z <= B.bmax.z) 
        return true;
    else return false;
}

bool segIntersect1d(const REAL & x1, const REAL & x2, const REAL & y1, const REAL & y2, REAL & r1, REAL & r2)
{
    /*

            y1        y2
            ------------
    -------------
    x1         x2

    */

    REAL A = x1 > y1 ? x1 : y1;
    REAL B = x2 < y2 ? x2 : y2;
    if (A > B) return false;
    else {
        r1 = A, r2 = B;
        return true;
    }
}

bool aabbIntersect(const aabb & A, const aabb & B, aabb & result)
{
    bool x = segIntersect1d(A.bmin.x, A.bmax.x, B.bmin.x, B.bmax.x, result.bmin.x, result.bmax.x);
    bool y = segIntersect1d(A.bmin.y, A.bmax.y, B.bmin.y, B.bmax.y, result.bmin.y, result.bmax.y);
    bool z = segIntersect1d(A.bmin.z, A.bmax.z, B.bmin.z, B.bmax.z, result.bmin.z, result.bmax.z);
    return x && y && z;
}

bool isIntersect(const triangle & A, const plane & B)
{
    VEC3 p1 = A.p[0] - B.p;
    VEC3 p2 = A.p[1] - B.p;
    VEC3 p3 = A.p[2] - B.p;
    bool sign1 = (dot(B.n, p1) > ZERO);
    bool sign2 = (dot(B.n, p2) > ZERO);
    bool sign3 = (dot(B.n, p3) > ZERO);
    if (sign1 && sign2 && sign3) return false;
    else if (!sign1 && !sign2 && !sign3) return false;
    else return true;
}

bool isIntersect(const sphere & A, const aabb & B)
{
    return distanceBetween(A.p, B) <= A.r ? true : false;
}

aabb computeBoundingBoxFromPointCloud(const Array<VEC3>& points)
{
    VEC3 _bmin = VEC3(REAL(REAL_MAX), REAL(REAL_MAX), REAL(REAL_MAX));
    VEC3 _bmax = VEC3(REAL(-REAL_MAX), REAL(-REAL_MAX), REAL(-REAL_MAX));

    for (int i = 0; i < points.size(); i++) {
        point p = points[i];
        if (_bmin.x > p.x) _bmin.x = p.x;
        if (_bmin.y > p.y) _bmin.y = p.y;
        if (_bmin.z > p.z) _bmin.z = p.z;
        if (_bmax.x < p.x) _bmax.x = p.x;
        if (_bmax.y < p.y) _bmax.y = p.y;
        if (_bmax.z < p.z) _bmax.z = p.z;
    }
    aabb bbox;
    bbox.bmin = _bmin;
    bbox.bmax = _bmax;
    return bbox;
}

int triangleAxisPlaneTest(const triangle & A, const char& axis, const REAL & value)
{
    int ax = axis - 'x';
    bool sign1 = (A.p[0].e[ax] > value);
    bool sign2 = (A.p[1].e[ax] > value);
    bool sign3 = (A.p[2].e[ax] > value);
    if (sign1 && sign2 && sign3) return 1;
    else if (!sign1 && !sign2 && !sign3) return -1;
    else return 0;
}

bool isIntersect(const triangle & A, const aabb & B)
{
    /* 
    when a triangle is "ON" the box, the function will return false,
    this is not intended so I added a small tolerance to avoid this
    from happening.
    */
    REAL tolerance = REAL(0.01);  /* expand box by 1% */
    VEC3 center = (B.bmax + B.bmin) * REAL(0.5);
    VEC3 half_size = (B.bmax - B.bmin) * REAL(0.5) * REAL(ONE + tolerance);

    REAL box_center[3];
    REAL box_halfsize[3];
    REAL triverts[3][3];

    box_center[0] = center.x; 
    box_center[1] = center.y; 
    box_center[2] = center.z;
    box_halfsize[0] = half_size.x; 
    box_halfsize[1] = half_size.y; 
    box_halfsize[2] = half_size.z;

    triverts[0][0] = A.p[0].x;
    triverts[0][1] = A.p[0].y;
    triverts[0][2] = A.p[0].z;
    triverts[1][0] = A.p[1].x;
    triverts[1][1] = A.p[1].y;
    triverts[1][2] = A.p[1].z;
    triverts[2][0] = A.p[2].x;
    triverts[2][1] = A.p[2].y;
    triverts[2][2] = A.p[2].z;

    if (_tri_box_overlap(box_center, box_halfsize, triverts)) return true;
    else return false;
}

bool isIntersect(const aabb & A, const aabb & B)
{
    /* x */
    REAL xLeft = A.bmin.x > B.bmin.x ? A.bmin.x : B.bmin.x;
    REAL xRight = A.bmax.x < B.bmax.x ? A.bmax.x : B.bmax.x;
    /* y */
    REAL yLeft = A.bmin.y > B.bmin.y ? A.bmin.y : B.bmin.y;
    REAL yRight = A.bmax.y < B.bmax.y ? A.bmax.y : B.bmax.y;
    /* z */
    REAL zLeft = A.bmin.z > B.bmin.z ? A.bmin.z : B.bmin.z;
    REAL zRight = A.bmax.z < B.bmax.z ? A.bmax.z : B.bmax.z;

    if (xLeft <= xRight && yLeft <= yRight && zLeft <= zRight)
        return true;
    else 
        return false;
}

bool isIntersect(const ray & A, const aabb & B)
{
    REAL t0 = ZERO, t1 = REAL_MAX, near_t, far_t;
    VEC3 invD = ONE / A.d;

    near_t = (B.bmin.x - A.o.x) * invD.x;
    far_t = (B.bmax.x - A.o.x) * invD.x;
    if (near_t > far_t) swap(near_t, far_t);
    t0 = t0 < near_t ? near_t : t0;
    t1 = t1 > far_t ? far_t : t1;
    if (t0 > t1) return false;

    near_t = (B.bmin.y - A.o.y) * invD.y;
    far_t = (B.bmax.y - A.o.y) * invD.y;
    if (near_t > far_t) swap(near_t, far_t);
    t0 = t0 < near_t ? near_t : t0;
    t1 = t1 > far_t ? far_t : t1;
    if (t0 > t1) return false;

    near_t = (B.bmin.z - A.o.z) * invD.z;
    far_t = (B.bmax.z - A.o.z) * invD.z;
    if (near_t > far_t) swap(near_t, far_t);
    t0 = t0 < near_t ? near_t : t0;
    t1 = t1 > far_t ? far_t : t1;
    if (t0 > t1) return false;

    return true;
}

bool isIntersect(const ray & A, const triangle & B)
{
    REAL u, v, t; /* ignored */
    return isIntersect(A, B, u, v, t);
}

bool isIntersect(const ray & A, const triangle & B, REAL& u, REAL& v, REAL& t)
{
    const REAL epsilon = REAL(1e-8);
    const REAL tolerance = REAL(1e-4);

    VEC3 e1 = B.p[1] - B.p[0];
    VEC3 e2 = B.p[2] - B.p[0];

    /* calculate plane's normal vector */
    VEC3 pvec = cross(A.d, e2);
    REAL det = dot(e1, pvec);
    if (Fabs(det) < epsilon) return false; /* ray is parallel to plane */

    REAL inv_det = ONE / det;
    VEC3 tvec = A.o - B.p[0];
    REAL i = dot(tvec, pvec) * inv_det;
    if (i < -tolerance || i > ONE + tolerance) return false;

    VEC3 qvec = cross(tvec, e1);
    REAL j = dot(A.d, qvec) * inv_det;
    if (j < -tolerance || i + j > ONE + tolerance) return false;
    t = dot(e2, qvec) * inv_det;
    if (t < -tolerance) return false;
    u = ONE - i - j;
    v = i;

    return true;
}

bool isIntersect(const triangle & A, const triangle & B)
{
    REAL v0[3], v1[3], v2[3];
    REAL u0[3], u1[3], u2[3];
    v0[0] = A.p[0].x; v0[1] = A.p[0].y; v0[2] = A.p[0].z;
    v1[0] = A.p[1].x; v1[1] = A.p[1].y; v1[2] = A.p[1].z;
    v2[0] = A.p[2].x; v2[1] = A.p[2].y; v2[2] = A.p[2].z;

    u0[0] = B.p[0].x; u0[1] = B.p[0].y; u0[2] = B.p[0].z;
    u1[0] = B.p[1].x; u1[1] = B.p[1].y; u1[2] = B.p[1].z;
    u2[0] = B.p[2].x; u2[1] = B.p[2].y; u2[2] = B.p[2].z;

    if (_tri_tri_intersect(v0, v1, v2, u0, u1, u2) > 0) return true;
    else return false;
}

bool isIntersect(const triangle & A, const triangle & B, int& isCoplanar, segment & s)
{
    REAL v0[3], v1[3], v2[3];
    REAL u0[3], u1[3], u2[3];
    v0[0] = A.p[0].x; v0[1] = A.p[0].y; v0[2] = A.p[0].z;
    v1[0] = A.p[1].x; v1[1] = A.p[1].y; v1[2] = A.p[1].z;
    v2[0] = A.p[2].x; v2[1] = A.p[2].y; v2[2] = A.p[2].z;

    u0[0] = B.p[0].x; u0[1] = B.p[0].y; u0[2] = B.p[0].z;
    u1[0] = B.p[1].x; u1[1] = B.p[1].y; u1[2] = B.p[1].z;
    u2[0] = B.p[2].x; u2[1] = B.p[2].y; u2[2] = B.p[2].z;

    REAL isect1[3], isect2[3];

    if (_tri_tri_intersect_with_isectline(v0, v1, v2, u0, u1, u2, 
        &isCoplanar, isect1, isect2) > 0) 
    {
        s.p[0].x = isect1[0]; s.p[0].y = isect1[1]; s.p[0].z = isect1[2];
        s.p[1].x = isect1[0]; s.p[1].y = isect1[1]; s.p[1].z = isect1[2];
        return true;
    }
    else return false;
}

bool isIntersect(const sphere & A, const triangle & B)
{
    if (A.r * A.r < squaredDistanceBetween(A.p, B)) 
        return false;
    else return true;
}

triangleEx::operator triangle() const
{
    triangle t;
    t.p[0] = this->p[0];
    t.p[1] = this->p[1];
    t.p[2] = this->p[2];
    return t;
}

ray::ray()
{
}

ray::ray(VEC3 o, VEC3 d)
{
    this->o = o;
    this->d = normalize(d);
}
