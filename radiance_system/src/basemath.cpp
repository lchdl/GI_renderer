#include "basemath.h"

REAL Cot(const REAL & x)
{
    return Tan(HALFPI - x);
}

VEC2 operator-(VEC2 a, VEC2 b)
{
    return VEC2(a.x - b.x, a.y - b.y);
}

VEC3 operator + (REAL a, VEC3 v) { return VEC3(a + v.x, a + v.y, a + v.z); }
VEC3 operator - (REAL a, VEC3 v) { return VEC3(a - v.x ,a - v.y, a - v.z); }
VEC3 operator * (REAL a, VEC3 v) { return VEC3(a * v.x, a * v.y, a * v.z); }
VEC3 operator / (REAL a, VEC3 v) { 
#ifdef USE_SAFE_DIVISION
    VEC3 v_s = v + REAL_MIN;
    return VEC3(a / v_s.x, a / v_s.y, a / v_s.z);
#else
    return VEC3(a / v.x, a / v.y, a / v.z);
#endif
}
VEC3 operator + (VEC3 v, REAL a) { return VEC3(v.x + a, v.y + a, v.z + a); }
VEC3 operator - (VEC3 v, REAL a) { return VEC3(v.x - a, v.y - a, v.z - a); }
VEC3 operator * (VEC3 v, REAL a) { return VEC3(v.x * a, v.y * a, v.z * a); }
VEC3 operator / (VEC3 v, REAL a) { 
#ifdef USE_SAFE_DIVISION
    REAL a_inv = ONE / (a + REAL_MIN);
#else
    REAL a_inv = ONE / a;
#endif
    return VEC3(v.x * a_inv, v.y * a_inv, v.z * a_inv);
}

VEC3 operator + (VEC3 a, VEC3 b) { return VEC3(a.x + b.x, a.y + b.y, a.z + b.z); }
VEC3 operator - (VEC3 a, VEC3 b) { return VEC3(a.x - b.x, a.y - b.y, a.z - b.z); }
VEC3 operator-(VEC3 v) { 
    return VEC3(-v.x, -v.y, -v.z); 
}

REAL dot(VEC3 a, VEC3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
VEC3 cross(VEC3 a, VEC3 b) {
    return VEC3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

REAL len(VEC3 a) { return Sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
REAL lensq(VEC3 a) { return (a.x * a.x + a.y * a.y + a.z * a.z); }
VEC3 normalize(VEC3 v) { 
    REAL l = len(v);
    if (Fabs(len(v)) > REAL_MIN) 
        return v / l;
    else {
#ifdef DEBUG
        printf("WARNING: Trying to normalize a zero vector.\n");
#endif
        return VEC3(ZERO, ZERO, ZERO);
    }
}

VEC3 operator*(VEC3 a, VEC3 b)
{
    return VEC3(a.x * b.x, a.y * b.y, a.z * b.z);
}

VEC4 operator+(REAL a, VEC4 v)
{
    return VEC4(a + v.s, a + v.x, a + v.y, a + v.z);
}

VEC4 operator-(REAL a, VEC4 v)
{
    return VEC4(a - v.s, a - v.x, a - v.y, a - v.z);
}

VEC4 operator*(REAL a, VEC4 v)
{
    return VEC4(a * v.s, a * v.x, a * v.y, a * v.z);
}

VEC4 operator/(REAL a, VEC4 v)
{
#ifdef USE_SAFE_DIVISION
    VEC4 v_s = v + REAL_MIN;
    return VEC4(a / v_s.s, a / v_s.x, a / v_s.y, a / v_s.z);
#else
    return VEC4(a / v.s, a / v.x, a / v.y, a / v.z);
#endif
}

VEC4 operator+(VEC4 v, REAL a)
{
    return VEC4(v.s + a, v.x + a, v.y + a, v.z + a);
}

VEC4 operator-(VEC4 v, REAL a)
{
    return VEC4(v.s - a, v.x - a, v.y - a, v.z - a);
}

VEC4 operator*(VEC4 v, REAL a)
{
    return VEC4(v.s * a, v.x * a, v.y * a, v.z * a);
}

VEC4 operator/(VEC4 v, REAL a)
{
#ifdef USE_SAFE_DIVISION
    REAL a_inv = ONE / (a + REAL_MIN);
#else
    REAL a_inv = ONE / a;
#endif
    return VEC4(v.s * a_inv, v.x * a_inv, v.y * a_inv, v.z * a_inv);
}

VEC4 operator/(VEC4 a, VEC4 b)
{
#ifdef USE_SAFE_DIVISION
    return VEC4(a.s / (b.s + REAL_MIN), a.x / (b.x + REAL_MIN), a.y / (b.y + REAL_MIN), a.z / (b.z + REAL_MIN));
#else
    return VEC4(a.s / b.s, a.x / b.x, a.y / b.y, a.z / b.z);
#endif
}

VEC4 VEC4::operator /= (const REAL& b)
{
#ifdef USE_SAFE_DIVISION
    REAL a;
    if (b == ZERO) a = b + REAL_MIN;
    else a = b;
    REAL inv_b = ONE / a;
    this->s *= inv_b;
    this->x *= inv_b;
    this->y *= inv_b;
    this->z *= inv_b;
    return (*this);
#else
    REAL inv_b = ONE / b;
    this->s *= inv_b;
    this->x *= inv_b;
    this->y *= inv_b;
    this->z *= inv_b;
    return (*this);
#endif
}

VEC4 operator+(VEC4 a, VEC4 b)
{
    return VEC4(a.s + b.s, a.x + b.x, a.y + b.y, a.z + b.z);
}

VEC4 operator+=(VEC4 a, VEC4 b)
{
    a.s += b.s;
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

VEC4 operator-(VEC4 a, VEC4 b)
{
    return VEC4(a.s - b.s, a.x - b.x, a.y - b.y, a.z - b.z);
}

VEC4 operator-=(VEC4 a, VEC4 b)
{
    a.s -= b.s;
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

VEC4 operator-(VEC4 v)
{
    return VEC4(-v.s, -v.x, -v.y, -v.z);
}

REAL dot(VEC4 a, VEC4 b)
{
    return a.s * b.s + a.x * b.x + a.y * b.y + a.z * b.z;
}

REAL len(VEC4 a)
{
    return Sqrt(lensq(a));
}

REAL lensq(VEC4 a)
{
    return a.s * a.s + a.x * a.x + a.y * a.y + a.z * a.z;
}

VEC4 normalize(VEC4 v)
{
    REAL l = len(v);
    if (Fabs(len(v)) > REAL_MIN)
        return v / l;
    else {
#ifdef DEBUG
        printf("WARNING: Trying to normalize a zero vector.\n");
#endif
        return VEC4(ZERO, ZERO, ZERO, ZERO);
    }
}

VEC4 operator*(VEC4 a, VEC4 b)
{
    return VEC4(a.s * b.s, a.x * b.x, a.y * b.y, a.z * b.z);
}

QUAT operator+(QUAT q1, QUAT q2) { return QUAT(q1.s + q2.s, q1.x + q2.x, q1.y + q2.y, q1.z + q2.z); }
QUAT operator-(QUAT q1, QUAT q2) { return QUAT(q1.s - q2.s, q1.x - q2.x, q1.y - q2.y, q1.z - q2.z); }

QUAT operator*(QUAT q1, QUAT q2)
{
    VEC3 v1 = VEC3(q1.x, q1.y, q1.z);
    VEC3 v2 = VEC3(q2.x, q2.y, q2.z);
    REAL s1 = q1.s;
    REAL s2 = q2.s;
    return QUAT(
        s1 * s2 - dot(v1, v2),
        s1 * v2 + s2 * v1 + cross(v1, v2)        
    );
}

QUAT operator/(QUAT q1, QUAT q2)
{
    return q1 * inv(q2);
}

QUAT conj(QUAT q)
{
    return QUAT(q.s, -q.x, -q.y, -q.z);
}

QUAT inv(QUAT q)
{
    REAL d = normsq(q);
#ifdef DEBUG
    if (Fabs(d) < REAL_MIN)
        printf("WARNING: quat inv(): division by zero.\n");
#endif
    return conj(q) / d;
}

REAL norm(QUAT q)
{
    return Sqrt( q.x * q.x + q.y * q.y + q.z * q.z + q.s * q.s);
}

REAL normsq(QUAT q) { return (q.x * q.x + q.y * q.y + q.z * q.z + q.s * q.s); }
QUAT normalize(QUAT q) { 
    REAL d = norm(q);
#ifdef DEBUG
    if (Fabs(d) < REAL_MIN)
        printf("WARNING: quat normalize(): division by zero.\n");
#endif

    return q / d;
}

REAL dot(QUAT q1, QUAT q2)
{
    return q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.s * q2.s;
}

QUAT slerp(QUAT q1, QUAT q2, REAL t)
{
    REAL theta = ArcCos(dot(q1, q2));
#ifdef DEBUG
    if (Fabs(Sin(theta)) < REAL_MIN)
        printf("WARNING: slerp(): division by zero.\n");
#endif
    QUAT q = (Sin((1 - t)*theta) * q1 + Sin(t*theta) * q2) / Sin(theta);
    return q;
}

MAT3x3 quatToMatrix(QUAT q)
{
    REAL t1, t2;
    MAT3x3 m;

    m.xx = ONE - TWO * (q.y * q.y + q.z * q.z);
    m.yy = ONE - TWO * (q.x * q.x + q.z * q.z);
    m.zz = ONE - TWO * (q.x * q.x + q.y * q.y);

    t1 = q.x * q.y;
    t2 = q.s * q.z;
    m.xy = TWO * (t1 - t2);
    m.yx = TWO * (t1 + t2);

    t1 = q.x * q.z;
    t2 = q.s * q.y;
    m.xz = TWO * (t1 + t2);
    m.zx = TWO * (t1 - t2);

    t1 = q.y * q.z;
    t2 = q.s * q.x;
    m.yz = TWO * (t1 - t2);
    m.zy = TWO * (t1 + t2);

    return m;
}

QUAT matrixToQuat(MAT3x3 m)
{
    REAL tr, s;
    int  i = 0;
    QUAT q;

    tr = m.xx + m.yy + m.zz;

    if (tr >= 0) {
        s = Sqrt(tr + ONE);
        q.s = HALF * s;
        s = HALF / s;
        q.x = (m.zy - m.yz) * s;
        q.y = (m.xz - m.zx) * s;
        q.z = (m.yx - m.xy) * s;
    }
    else {
        if (m.yy > m.xx) i = 4;
        if (m.zz > m.e[i]) i = 8;
        switch (i) {
        case 0:
            s = Sqrt((m.xx - (m.yy + m.zz)) + ONE);
            q.x = HALF * s;
            s = HALF / s;
            q.y = (m.xy + m.yx) * s;
            q.z = (m.zx + m.xz) * s;
            q.s = (m.zy - m.yz) * s;
            break;
        case 4:
            s = Sqrt((m.yy - (m.zz + m.xx)) + ONE);
            q.y = HALF * s;
            s = HALF / s;
            q.z = (m.yz + m.zy) * s;
            q.x = (m.xy + m.yx) * s;
            q.s = (m.xz - m.zx) * s;
            break;
        case 8:
            s = Sqrt((m.zz - (m.xx + m.yy)) + ONE);
            q.z = HALF * s;
            s = HALF / s;
            q.x = (m.zx + m.xz) * s;
            q.y = (m.yz + m.zy) * s;
            q.s = (m.yx - m.xy) * s;
            break;
        }
    }

    return q;
}

QUAT vecMulQuat(VEC3 v, QUAT q)
{
    QUAT t(ZERO, v); /* expand vector to quaternion */
    return t * q;
}

QUAT eulerToQuat(REAL yaw, REAL pitch, REAL roll)
{
    /* https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles */

    REAL cy = Cos(yaw * HALF);
    REAL sy = Sin(yaw * HALF);
    REAL cp = Cos(pitch * HALF);
    REAL sp = Sin(pitch * HALF);
    REAL cr = Cos(roll * HALF);
    REAL sr = Sin(roll * HALF);

    return QUAT(
        cr * cp * cy + sr * sp * sy,
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy);
}

MAT3x3 operator*(MAT3x3 a, MAT3x3 b)
{
    return MAT3x3(
        a.xx * b.xx + a.xy * b.yx + a.xz * b.zx, a.xx * b.xy + a.xy * b.yy + a.xz * b.zy, a.xx * b.xz + a.xy * b.yz + a.xz * b.zz,
        a.yx * b.xx + a.yy * b.yx + a.yz * b.zx, a.yx * b.xy + a.yy * b.yy + a.yz * b.zy, a.yx * b.xz + a.yy * b.yz + a.yz * b.zz,
        a.zx * b.xx + a.zy * b.yx + a.zz * b.zx, a.zx * b.xy + a.zy * b.yy + a.zz * b.zy, a.zx * b.xz + a.zy * b.yz + a.zz * b.zz
    );
}

VEC3 operator*(MAT3x3 a, VEC3 b)
{
    return VEC3(
        a.xx * b.x + a.xy * b.y + a.xz * b.z,
        a.yx * b.x + a.yy * b.y + a.yz * b.z,
        a.zx * b.x + a.zy * b.y + a.zz * b.z
    );
}

MAT3x3 transpose(MAT3x3 a)
{
    return MAT3x3(
        a.xx, a.yx, a.zx,
        a.xy, a.yy, a.zy,
        a.xz, a.yz, a.zz
    );
}

MAT3x3 inv(MAT3x3 a)
{
	MAT3x3 aInv;
	REAL det = a.xx * (a.yy * a.zz - a.zy * a.yz) -
		a.xy * (a.yx * a.zz - a.yz * a.zx) +
		a.xz * (a.yx * a.zy - a.yy * a.zx);
	REAL detInv = ONE / det;
#ifdef DEBUG
    if (Fabs(det) < REAL_MIN)
        printf("WARNING: Matrix cannot be inversed (det=0).\n");
#endif
	aInv.xx = (a.yy * a.zz - a.zy * a.yz) * detInv;
	aInv.xy = (a.xz * a.zy - a.xy * a.zz) * detInv;
	aInv.xz = (a.xy * a.yz - a.xz * a.yy) * detInv;
	aInv.yx = (a.yz * a.zx - a.yx * a.zz) * detInv;
	aInv.yy = (a.xx * a.zz - a.xz * a.zx) * detInv;
	aInv.yz = (a.yx * a.xz - a.xx * a.yz) * detInv;
	aInv.zx = (a.yx * a.zy - a.zx * a.yy) * detInv;
	aInv.zy = (a.zx * a.xy - a.xx * a.zy) * detInv;
	aInv.zz = (a.xx * a.yy - a.yx * a.xy) * detInv;
	return aInv;
}

MAT4x4 transpose(MAT4x4 a)
{
    return MAT4x4(
        a.xx, a.yx, a.zx, a.sx,
        a.xy, a.yy, a.zy, a.sy,
        a.xz, a.yz, a.zz, a.sz,
        a.xs, a.ys, a.zs, a.ss
    );
}

MAT3x3 operator*(REAL a, MAT3x3 m)
{
    return MAT3x3(
        a * m.xx, a * m.xy, a * m.xz,
        a * m.yx, a * m.yy, a * m.yz,
        a * m.zx, a * m.zy, a * m.zz
    );
}

MAT3x3 operator*(MAT3x3 m, REAL a)
{
    return MAT3x3(
        a * m.xx, a * m.xy, a * m.xz,
        a * m.yx, a * m.yy, a * m.yz,
        a * m.zx, a * m.zy, a * m.zz
    );
}

VEC3 VecProject(VEC3 v, VEC3 n)
{
    return dot(v, n) * n;
}

VEC3 VecTangent(VEC3 v, VEC3 n)
{
    return v - VecProject(v, n);
}


/* https://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi */
/* wrap x -> [0,max) */
REAL wrap(REAL x, REAL max)
{
    /* integer math: `(max + x % max) % max` */
    return Fmod(max + Fmod(x, max), max);
}

/* wrap x -> [min,max) */
REAL wrap(REAL x, REAL min, REAL max)
{
    return min + wrap(x - min, max - min);
}

REAL degToRad(const REAL & deg)
{
    return deg * PI / REAL(180.0);
}

REAL radToDeg(const REAL & rad)
{
    return rad / PI * REAL(180.0);
}

QUAT operator*(QUAT q, REAL a) { return QUAT(q.s * a, q.x * a, q.y * a, q.z * a); }
QUAT operator/(QUAT q, REAL a) { 
#ifdef USE_SAFE_DIVISION
    REAL a_s = a + REAL_MIN;
    return QUAT(q.s / a_s, q.x / a_s, q.y / a_s, q.z / a_s);
#else
    return QUAT(q.s / a, q.x / a, q.y / a, q.z / a);
#endif
}
QUAT operator*(REAL a, QUAT q) { return QUAT(a * q.s, a * q.x, a * q.y, a * q.z); }
QUAT operator/(REAL a, QUAT q) { return a * inv(q); }

INT3 operator+(INT3 a, INT3 b) { return INT3(a.x + b.x, a.y + b.y, a.z + b.z); }
INT3 operator-(INT3 a, INT3 b) { return INT3(a.x - b.x, a.y - b.y, a.z - b.z); }
INT3 operator*(INT3 a, INT3 b) { return INT3(a.x * b.x, a.y * b.y, a.z * b.z); }
INT3 operator/(INT3 a, INT3 b) { return INT3(a.x / b.x, a.y / b.y, a.z / b.z); }

bool operator==(VEC3 a, VEC3 b)
{
    if (len(a - b) <= 0.02) return true;
    else return false;
}

VEC4 VEC4::operator+=(const VEC4 & a)
{
    this->s += a.s;
    this->x += a.x;
    this->y += a.y;
    this->z += a.z;
    return *this;
}

VEC4 VEC4::operator-=(const VEC4 & a)
{
    this->s -= a.s;
    this->x -= a.x;
    this->y -= a.y;
    this->z -= a.z;
    return *this;
}

VEC3 VEC3::operator+=(const VEC3 & v)
{
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return (*this);
}

VEC3 VEC3::operator-=(const VEC3 & v)
{
    this->x -= v.x;
    this->y -= v.y;
    this->z -= v.z;
    return (*this);
}

VEC3 VEC3::operator*=(const VEC3 & v)
{
    this->x *= v.x;
    this->y *= v.y;
    this->z *= v.z;
    return (*this);
}

VEC3 VEC3::operator/=(const VEC3 & v)
{
#ifdef USE_SAFE_DIVISION
    VEC3 div = VEC3(
        v.x == ZERO ? REAL_MIN : v.x,
        v.y == ZERO ? REAL_MIN : v.y,
        v.z == ZERO ? REAL_MIN : v.z);
    this->x /= div.x;
    this->y /= div.y;
    this->z /= div.z;
    return (*this);
#else
    this->x /= v.x;
    this->y /= v.y;
    this->z /= v.z;
    return (*this);
#endif
}

VEC3 VEC3::operator*=(const REAL & a)
{
    this->x *= a;
    this->y *= a;
    this->z *= a;
    return (*this);
}

VEC3 VEC3::operator/=(const REAL & a)
{
#ifdef USE_SAFE_DIVISION
    REAL div = (a == ZERO ? REAL_MIN : a);
#else
    REAL div = a;
#endif
    REAL inv_div = ONE / div;
    this->x *= inv_div;
    this->y *= inv_div;
    this->z *= inv_div;
    return (*this);
}
