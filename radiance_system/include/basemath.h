#pragma once

#include "global.h"
#include <math.h>
#include <float.h>
#include <stdio.h>

/* define division by zero error */
/* Current supported platforms: WINDOWS or LINUX-based systems. */
#if defined(_MSC_VER) /* if compile this program on MSVC/VS IDE */
#define RAISE_DIVISION_BY_ZERO_ERROR __debugbreak(); /* then use builtin function */
#elif defined(__linux__)
#include <signal.h>
#define RAISE_DIVISION_BY_ZERO_ERROR raise(SIGTRAP); /* else we use raise(...) */
#else
#error "Cannot define RAISE_DIVISION_BY_ZERO_ERROR macro due to unknown platform or compiler."
#endif

#ifdef USE_DOUBLE_PRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

/* Define unified macros for single/double precision math functions */
#ifdef USE_DOUBLE_PRECISION
#define Exp    exp
#define Sqrt   sqrt
#define Sin    sin
#define Cos    cos
#define ArcCos acos
#define Pow    pow
#define Fmod   fmod
#define Modf   modf
#define Fabs   fabs
#define Round  round
#define Log    log
#define Tan    tan
#define Arctan atan
#define Arccos acos
#define Arcsin asin
#else /* USE_SINGLE_PRECISION */
#define Exp    expf
#define Sqrt   sqrtf
#define Sin    sinf
#define Cos    cosf
#define ArcCos acosf
#define Pow    powf
#define Fmod   fmodf
#define Modf   modff
#define Fabs   fabsf
#define Round  roundf
#define Log    logf
#define Tan    tanf
#define Arctan atanf
#define Arccos acosf
#define Arcsin asinf
#endif
#define Ln Log

#ifdef USE_DOUBLE_PRECISION
#define REAL_MIN DBL_MIN
#define REAL_MAX DBL_MAX
#else
#define REAL_MIN FLT_MIN
#define REAL_MAX FLT_MAX
#endif

#define ZERO    ((REAL)(0.0))  /* cast constants to proper type (double or float) */
#define QUARTER ((REAL)(0.25))
#define HALF    ((REAL)(0.5))
#define ONE     ((REAL)(1.0))
#define TWO     ((REAL)(2.0))

#define HALFPI  ((REAL)(1.57079632679))
#define PI      ((REAL)(3.14159265359))
#define TWOPI   ((REAL)(6.28318530718))
#define PI_LP   ((REAL)(3.14)) /* low-precision PI */
#define PI_INV  ((REAL)(0.31830988618))

REAL Cot(const REAL& x); /* cotangent */

#define MIN(a,b) (( (a) < (b) ) ? (a) : (b) )
#define MAX(a,b) (( (a) > (b) ) ? (a) : (b) )

/* vector with two elements */
struct VEC2 {
    union {
        REAL e[2];      /* packed components */
        struct {
            REAL x, y;  /* unpacked components */
        };
    };

    VEC2() { x = ZERO; y = ZERO; }
    VEC2(REAL x, REAL y) { this->x = x; this->y = y; }
    static VEC2 ones() { return VEC2(ONE, ONE); }
    static VEC2 zeros() { return VEC2(ZERO, ZERO); }
};

/* vector with three elements */
struct VEC3 {
    union {
        REAL e[3];
        struct {
            REAL x, y, z;
        };
    };

    VEC3() { x = ZERO; y = ZERO; z = ZERO; }
    VEC3(REAL x, REAL y, REAL z) { this->x = x; this->y = y; this->z = z; }

    static VEC3 ones() { return VEC3(ONE, ONE, ONE); }
    static VEC3 zeros() { return VEC3(ZERO, ZERO, ZERO); }
    static VEC3 all(REAL a) { return VEC3(a, a, a); }
    static VEC3 min() { return VEC3(-REAL_MAX, -REAL_MAX, -REAL_MAX); }
    static VEC3 max() { return VEC3(REAL_MAX, REAL_MAX, REAL_MAX); }

    VEC3 operator += (const VEC3& v);
    VEC3 operator -= (const VEC3& v);
    VEC3 operator *= (const VEC3& v);
    VEC3 operator /= (const VEC3& v);
    VEC3 operator *= (const REAL& a);
    VEC3 operator /= (const REAL& a);

};

struct VEC4 {
    union {
        REAL e[4];
        struct {
            REAL s, x, y, z;
        };
    };
    VEC4() { s = ZERO; x = ZERO; y = ZERO; z = ZERO; }
    VEC4(REAL s, REAL x, REAL y, REAL z) { this->s = s; this->x = x; this->y = y; this->z = z; }
    VEC4(VEC3 v, REAL a) { this->s = v.x; this->x = v.y; this->y = v.z; this->z = a; }
    VEC4(REAL a, VEC3 v) { this->s = a; this->x = v.x; this->y = v.y; this->z = v.z; }
    VEC3 xyz() const { return VEC3(x, y, z); }

    static VEC4 ones() { return VEC4(ONE, ONE, ONE, ONE); }
    static VEC4 zeros() { return VEC4(ZERO, ZERO, ZERO, ZERO); }
    static VEC4 all(REAL a) { return VEC4(a, a, a, a); }
    static VEC4 min() { return VEC4(-REAL_MAX, -REAL_MAX, -REAL_MAX, -REAL_MAX); }
    static VEC4 max() { return VEC4(REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX); }

    VEC4 operator += (const VEC4& a);
    VEC4 operator -= (const VEC4& a);
    VEC4 operator /= (const REAL& b);

};

/* Quaternion. q = s + xi + yj + zk */
struct QUAT {
    union {
        REAL e[4];
        struct {
            REAL s, x, y, z;
        };
    };

    QUAT() { x = ZERO; y = ZERO; z = ZERO; s = ZERO; }
    QUAT(REAL s, REAL x, REAL y, REAL z) { this->s = s; this->x = x; this->y = y; this->z = z; }
    QUAT(REAL s, VEC3 v) { this->s = s; this->x = v.x; this->y = v.y; this->z = v.z; }

    static QUAT identity() { return QUAT(ONE, VEC3::zeros()); }
    static QUAT rotX(REAL angle) { return QUAT(Cos(angle / TWO), Sin(angle / TWO), ZERO, ZERO); }
    static QUAT rotY(REAL angle) { return QUAT(Cos(angle / TWO), ZERO, Sin(angle / TWO), ZERO); }
    static QUAT rotZ(REAL angle) { return QUAT(Cos(angle / TWO), ZERO, ZERO, Sin(angle / TWO)); }

};

/* a 3x3 matrix */
struct MAT3x3 {
    union {
        REAL e[9];
        struct {
            REAL xx, xy, xz;
            REAL yx, yy, yz;
            REAL zx, zy, zz;
        };
    };

    MAT3x3() { for (int i = 0; i < 9; i++) e[i] = ZERO; }
    MAT3x3(
        REAL xx, REAL xy, REAL xz,
        REAL yx, REAL yy, REAL yz,
        REAL zx, REAL zy, REAL zz) {
        e[0] = xx, e[1] = xy, e[2] = xz;
        e[3] = yx, e[4] = yy, e[5] = yz;
        e[6] = zx, e[7] = zy, e[8] = zz;
    }
    static MAT3x3 identity() {
        return MAT3x3(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1);
    }
    static MAT3x3 diag(const REAL& xx, const REAL& yy, const REAL& zz) {
        return MAT3x3(
            xx, 0, 0,
            0, yy, 0,
            0, 0, zz);
    }
};
typedef MAT3x3 MAT3;

struct MAT4x4 {
    union {
        REAL e[16];
        struct {
            REAL xx, xy, xz, xs;
            REAL yx, yy, yz, ys;
            REAL zx, zy, zz, zs;
            REAL sx, sy, sz, ss;
        };
    };

    MAT4x4() { for (int i = 0; i < 16; i++) e[i] = ZERO; }
    MAT4x4(
        REAL xx, REAL xy, REAL xz, REAL xs,
        REAL yx, REAL yy, REAL yz, REAL ys,
        REAL zx, REAL zy, REAL zz, REAL zs,
        REAL sx, REAL sy, REAL sz, REAL ss) {
        e[0] = xx, e[1] = xy, e[2] = xz; e[3] = xs;
        e[4] = yx, e[5] = yy, e[6] = yz; e[7] = ys;
        e[8] = zx, e[9] = zy, e[10] = zz; e[11] = zs;
        e[12] = sx, e[13] = sy, e[14] = sz; e[15] = ss;
    }
    static MAT4x4 identity() {
        return MAT4x4(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1);
    }
};
typedef MAT4x4 MAT4;

struct INT2 {
    union {
        int e[2];
        struct {
            int x, y;
        };
    };
    INT2() { x = 0; y = 0; }
    INT2(int x, int y) { this->x = x; this->y = y; }
};

struct INT3 {
    union {
        int e[3];
        struct {
            int x, y, z;
        };
    };
    INT3() { x = 0; y = 0; z = 0; }
    INT3(int x, int y, int z) { this->x = x; this->y = y; this->z = z; }
};

/* basic mathematical operations for vectors, quaternions and matrices */
VEC2 operator - (VEC2 a, VEC2 b);

VEC3 operator + (REAL a, VEC3 v);
VEC3 operator - (REAL a, VEC3 v);
VEC3 operator * (REAL a, VEC3 v);
VEC3 operator / (REAL a, VEC3 v);
VEC3 operator + (VEC3 v, REAL a);
VEC3 operator - (VEC3 v, REAL a);
VEC3 operator * (VEC3 v, REAL a);
VEC3 operator / (VEC3 v, REAL a);

VEC3 operator + (VEC3 a, VEC3 b);
VEC3 operator - (VEC3 a, VEC3 b);
VEC3 operator - (VEC3 v); /* negate */
REAL dot(VEC3 a, VEC3 b);
VEC3 cross(VEC3 a, VEC3 b);
REAL len(VEC3 a);
REAL lensq(VEC3 a);
VEC3 normalize(VEC3 v);

VEC3 operator * (VEC3 a, VEC3 b);

VEC4 operator + (REAL a, VEC4 v);
VEC4 operator - (REAL a, VEC4 v);
VEC4 operator * (REAL a, VEC4 v);
VEC4 operator / (REAL a, VEC4 v);
VEC4 operator + (VEC4 v, REAL a);
VEC4 operator - (VEC4 v, REAL a);
VEC4 operator * (VEC4 v, REAL a);
VEC4 operator / (VEC4 v, REAL a);
VEC4 operator / (VEC4 a, VEC4 b);

VEC4 operator + (VEC4 a, VEC4 b);
VEC4 operator - (VEC4 a, VEC4 b);
VEC4 operator - (VEC4 v); /* negate */

REAL dot(VEC4 a, VEC4 b);
REAL len(VEC4 a);
REAL lensq(VEC4 a);
VEC4 normalize(VEC4 v);

VEC4 operator * (VEC4 a, VEC4 b);

QUAT operator + (QUAT q1, QUAT q2);
QUAT operator - (QUAT q1, QUAT q2);
QUAT operator * (QUAT q1, QUAT q2);    /* Hamilton product of two quaternions */
QUAT operator * (QUAT q, REAL a);
QUAT operator * (REAL a, QUAT q);
QUAT operator / (QUAT q1, QUAT q2);
QUAT operator / (QUAT q, REAL a);
QUAT operator / (REAL a, QUAT q);

INT3 operator + (INT3 a, INT3 b);
INT3 operator - (INT3 a, INT3 b);
INT3 operator * (INT3 a, INT3 b);
INT3 operator / (INT3 a, INT3 b);

/* for debugging ONLY */
bool operator == (VEC3 a, VEC3 b);

QUAT conj(QUAT q);                           /* conjugate */
QUAT inv(QUAT q);                            /* inverse */
REAL norm(QUAT q);                           /* returns the norm ||q|| of a quaternion */
REAL normsq(QUAT q);                         /* returns squared norm  */
QUAT normalize(QUAT q);                      /* normalizing a quaternion */
REAL dot(QUAT q1, QUAT q2);
QUAT slerp(QUAT q1, QUAT q2, REAL t);  /* spherical linear interpolation */

MAT3x3 quatToMatrix(QUAT q);
QUAT matrixToQuat(MAT3x3 m);
QUAT vecMulQuat(VEC3 v, QUAT q);

QUAT eulerToQuat(REAL yaw, REAL pitch, REAL roll);

MAT3x3 operator * (MAT3x3 a, MAT3x3 b);
VEC3 operator * (MAT3x3 a, VEC3 b);
MAT3x3 transpose(MAT3x3 a);
MAT3x3 inv(MAT3x3 a);

MAT4x4 transpose(MAT4x4 a);

MAT3x3 operator * (REAL a, MAT3x3 m);
MAT3x3 operator * (MAT3x3 m, REAL a);

VEC3 VecProject(VEC3 v, VEC3 n);
VEC3 VecTangent(VEC3 v, VEC3 n);

template <typename T>
T max2(const T& a, const T& b)
{
    return a > b ? a : b;
}
template <typename T>
T min2(const T& a, const T& b)
{
    return a < b ? a : b;
}
template <typename T>
T max3(const T& a, const T& b, const T& c)
{
    T m = a > b ? a : b;
    if (m < c) m = c;
    return m;
}
template <typename T>
int argmax3(const T& a, const T& b, const T& c)
{
    int arg = 0;
    if (a > b) {
        if (a > c) arg = 0;
        else arg = 2;
    }
    else {
        if (b > c) arg = 1;
        else arg = 2;
    }
    return arg;
}
template <typename T>
T min3(const T& a, const T& b, const T& c)
{
    T m = a < b ? a : b;
    if (m > c) m = c;
    return m;
}
template <typename T>
int argmin3(const T& a, const T& b, const T& c)
{
    int arg = 0;
    if (a < b) {
        if (a < c) arg = 0;
        else arg = 2;
    }
    else {
        if (b < c) arg = 1;
        else arg = 2;
    }
    return arg;
}
template <typename T>
int argmin(const T* a, const int n, T* pmin = NULL)
{
    int arg = 0;
    T m = a[0];
    for (int i = 1; i < n; i++) {
        if (m > a[i]) {
            m = a[i];
            arg = i;
        }
    }
    if (pmin != NULL)
        *pmin = m;
    return arg;
}
template <typename T>
int rangedArgmin(const T* a, const int start, const int end, T* pmin = NULL)
{
    int arg = start;
    T m = a[start];
    for (int i = start; i < end; i++) {
        if (m > a[i]) {
            m = a[i];
            arg = i;
        }
    }
    if (pmin != NULL)
        *pmin = m;
    return arg;
}
template <typename T>
T rangedMin(const T* a, const int start, const int end)
{
    T m = a[start];
    for (int i = start; i < end; i++) {
        if (m > a[i]) {
            m = a[i];
        }
    }
    return m;
}

/* https://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi */
/* wrap x -> [0,max) */
REAL wrap(REAL x, REAL max);
/* wrap x -> [min,max) */
REAL wrap(REAL x, REAL min, REAL max);

REAL degToRad(const REAL& deg);
REAL radToDeg(const REAL& rad);

