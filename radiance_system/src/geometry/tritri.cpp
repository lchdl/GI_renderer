#include "geometry/tritri.h"

int _coplanar_tri_tri(REAL N[3], REAL V0[3], REAL V1[3], REAL V2[3], REAL U0[3], REAL U1[3], REAL U2[3])
{
    REAL A[3];
    short i0, i1;
    /* first project onto an axis-aligned plane, that maximizes the area */
    /* of the triangles, compute indices: i0,i1. */
    A[0] = Fabs(N[0]);
    A[1] = Fabs(N[1]);
    A[2] = Fabs(N[2]);
    if (A[0] > A[1])
    {
        if (A[0] > A[2])
        {
            i0 = 1;      /* A[0] is greatest */
            i1 = 2;
        }
        else
        {
            i0 = 0;      /* A[2] is greatest */
            i1 = 1;
        }
    }
    else   /* A[0]<=A[1] */
    {
        if (A[2] > A[1])
        {
            i0 = 0;      /* A[2] is greatest */
            i1 = 1;
        }
        else
        {
            i0 = 0;      /* A[1] is greatest */
            i1 = 2;
        }
    }

    /* test all edges of triangle 1 against the edges of triangle 2 */
    _TRITRI_EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2);
    _TRITRI_EDGE_AGAINST_TRI_EDGES(V1, V2, U0, U1, U2);
    _TRITRI_EDGE_AGAINST_TRI_EDGES(V2, V0, U0, U1, U2);

    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    _TRITRI_POINT_IN_TRI(V0, U0, U1, U2);
    _TRITRI_POINT_IN_TRI(U0, V0, V1, V2);

    return 0;
}

int _tri_tri_intersect(REAL V0[3], REAL V1[3], REAL V2[3], REAL U0[3], REAL U1[3], REAL U2[3])
{
    REAL E1[3], E2[3];
    REAL N1[3], N2[3], d1, d2;
    REAL du0, du1, du2, dv0, dv1, dv2;
    REAL D[3];
    REAL isect1[2], isect2[2];
    REAL du0du1, du0du2, dv0dv1, dv0dv2;
    short index;
    REAL vp0, vp1, vp2;
    REAL up0, up1, up2;
    REAL bb, cc, max;
    REAL a, b, c, x0, x1;
    REAL d, e, f, y0, y1;
    REAL xx, yy, xxyy, tmp;

    /* compute plane equation of triangle(V0,V1,V2) */
    _TRITRI_SUB(E1, V1, V0);
    _TRITRI_SUB(E2, V2, V0);
    _TRITRI_CROSS(N1, E1, E2);
    d1 = -_TRITRI_DOT(N1, V0);
    /* plane equation 1: N1.X+d1=0 */

    /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
    du0 = _TRITRI_DOT(N1, U0) + d1;
    du1 = _TRITRI_DOT(N1, U1) + d1;
    du2 = _TRITRI_DOT(N1, U2) + d1;

    /* coplanarity robustness check */
#if _USE_EPSILON_TEST==TRUE
    if (Fabs(du0) < _TRITRI_EPSILON) du0 = ZERO;
    if (Fabs(du1) < _TRITRI_EPSILON) du1 = ZERO;
    if (Fabs(du2) < _TRITRI_EPSILON) du2 = ZERO;
#endif
    du0du1 = du0 * du1;
    du0du2 = du0 * du2;

    if (du0du1 > ZERO && du0du2 > ZERO) /* same sign on all of them + not equal 0 ? */
        return 0;                    /* no intersection occurs */

                                     /* compute plane of triangle (U0,U1,U2) */
    _TRITRI_SUB(E1, U1, U0);
    _TRITRI_SUB(E2, U2, U0);
    _TRITRI_CROSS(N2, E1, E2);
    d2 = -_TRITRI_DOT(N2, U0);
    /* plane equation 2: N2.X+d2=0 */

    /* put V0,V1,V2 into plane equation 2 */
    dv0 = _TRITRI_DOT(N2, V0) + d2;
    dv1 = _TRITRI_DOT(N2, V1) + d2;
    dv2 = _TRITRI_DOT(N2, V2) + d2;

#if _USE_EPSILON_TEST==TRUE
    if (Fabs(dv0) < _TRITRI_EPSILON) dv0 = ZERO;
    if (Fabs(dv1) < _TRITRI_EPSILON) dv1 = ZERO;
    if (Fabs(dv2) < _TRITRI_EPSILON) dv2 = ZERO;
#endif

    dv0dv1 = dv0 * dv1;
    dv0dv2 = dv0 * dv2;

    if (dv0dv1 > ZERO && dv0dv2 > ZERO) /* same sign on all of them + not equal 0 ? */
        return 0;                    /* no intersection occurs */

                                     /* compute direction of intersection line */
    _TRITRI_CROSS(D, N1, N2);

    /* compute and index to the largest component of D */
    max = (REAL)Fabs(D[0]);
    index = 0;
    bb = (REAL)Fabs(D[1]);
    cc = (REAL)Fabs(D[2]);
    if (bb > max) max = bb, index = 1;
    if (cc > max) max = cc, index = 2;

    /* this is the simplified projection onto L*/
    vp0 = V0[index];
    vp1 = V1[index];
    vp2 = V2[index];

    up0 = U0[index];
    up1 = U1[index];
    up2 = U2[index];

    /* compute interval for triangle 1 */
    _TRITRI_NEWCOMPUTE_INTERVALS(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, a, b, c, x0, x1);

    /* compute interval for triangle 2 */
    _TRITRI_NEWCOMPUTE_INTERVALS(up0, up1, up2, du0, du1, du2, du0du1, du0du2, d, e, f, y0, y1);

    xx = x0 * x1;
    yy = y0 * y1;
    xxyy = xx * yy;

    tmp = a * xxyy;
    isect1[0] = tmp + b * x1*yy;
    isect1[1] = tmp + c * x0*yy;

    tmp = d * xxyy;
    isect2[0] = tmp + e * xx*y1;
    isect2[1] = tmp + f * xx*y0;

    _TRITRI_SORT(isect1[0], isect1[1]);
    _TRITRI_SORT(isect2[0], isect2[1]);

    if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;
    return 1;
}

void _isect2(REAL VTX0[3], REAL VTX1[3], REAL VTX2[3], REAL VV0, REAL VV1, REAL VV2, REAL D0, REAL D1, REAL D2, REAL * isect0, REAL * isect1, REAL isectpoint0[3], REAL isectpoint1[3])
{
    REAL tmp = D0 / (D0 - D1);
    REAL diff[3];
    *isect0 = VV0 + (VV1 - VV0)*tmp;
    _TRITRI_SUB(diff, VTX1, VTX0);
    _TRITRI_MULT(diff, diff, tmp);
    _TRITRI_ADD(isectpoint0, diff, VTX0);
    tmp = D0 / (D0 - D2);
    *isect1 = VV0 + (VV2 - VV0)*tmp;
    _TRITRI_SUB(diff, VTX2, VTX0);
    _TRITRI_MULT(diff, diff, tmp);
    _TRITRI_ADD(isectpoint1, VTX0, diff);
}

int _compute_intervals_isectline(REAL VERT0[3], REAL VERT1[3], REAL VERT2[3], REAL VV0, REAL VV1, REAL VV2, REAL D0, REAL D1, REAL D2, REAL D0D1, REAL D0D2, REAL * isect0, REAL * isect1, REAL isectpoint0[3], REAL isectpoint1[3])
{
    if (D0D1 > ZERO)
    {
        /* here we know that D0D2<=ZERO */
        /* that is D0, D1 are on the same side, D2 on the other or on the plane */
        _isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, isect0, isect1, isectpoint0, isectpoint1);
    }
    else if (D0D2 > ZERO)
    {
        /* here we know that d0d1<=ZERO */
        _isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, isect0, isect1, isectpoint0, isectpoint1);
    }
    else if (D1*D2 > ZERO || D0 != ZERO)
    {
        /* here we know that d0d1<=ZERO or that D0!=ZERO */
        _isect2(VERT0, VERT1, VERT2, VV0, VV1, VV2, D0, D1, D2, isect0, isect1, isectpoint0, isectpoint1);
    }
    else if (D1 != ZERO)
    {
        _isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, isect0, isect1, isectpoint0, isectpoint1);
    }
    else if (D2 != ZERO)
    {
        _isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, isect0, isect1, isectpoint0, isectpoint1);
    }
    else
    {
        /* triangles are coplanar */
        return 1;
    }
    return 0;
}

int _tri_tri_intersect_with_isectline(REAL V0[3], REAL V1[3], REAL V2[3], REAL U0[3], REAL U1[3], REAL U2[3], int * coplanar, REAL isectpt1[3], REAL isectpt2[3])
{
    REAL E1[3], E2[3];
    REAL N1[3], N2[3], d1, d2;
    REAL du0, du1, du2, dv0, dv1, dv2;
    REAL D[3];
    REAL isect1[2], isect2[2];
    REAL isectpointA1[3], isectpointA2[3];
    REAL isectpointB1[3], isectpointB2[3];
    REAL du0du1, du0du2, dv0dv1, dv0dv2;
    short index;
    REAL vp0, vp1, vp2;
    REAL up0, up1, up2;
    REAL b, c, max;
    int smallest1, smallest2;

    /* compute plane equation of triangle(V0,V1,V2) */
    _TRITRI_SUB(E1, V1, V0);
    _TRITRI_SUB(E2, V2, V0);
    _TRITRI_CROSS(N1, E1, E2);
    d1 = -_TRITRI_DOT(N1, V0);
    /* plane equation 1: N1.X+d1=0 */

    /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
    du0 = _TRITRI_DOT(N1, U0) + d1;
    du1 = _TRITRI_DOT(N1, U1) + d1;
    du2 = _TRITRI_DOT(N1, U2) + d1;

    /* coplanarity robustness check */
#if _USE_EPSILON_TEST==TRUE
    if (Fabs(du0) < _TRITRI_EPSILON) du0 = ZERO;
    if (Fabs(du1) < _TRITRI_EPSILON) du1 = ZERO;
    if (Fabs(du2) < _TRITRI_EPSILON) du2 = ZERO;
#endif
    du0du1 = du0 * du1;
    du0du2 = du0 * du2;

    if (du0du1 > ZERO && du0du2 > ZERO) /* same sign on all of them + not equal 0 ? */
        return 0;                    /* no intersection occurs */

                                     /* compute plane of triangle (U0,U1,U2) */
    _TRITRI_SUB(E1, U1, U0);
    _TRITRI_SUB(E2, U2, U0);
    _TRITRI_CROSS(N2, E1, E2);
    d2 = -_TRITRI_DOT(N2, U0);
    /* plane equation 2: N2.X+d2=0 */

    /* put V0,V1,V2 into plane equation 2 */
    dv0 = _TRITRI_DOT(N2, V0) + d2;
    dv1 = _TRITRI_DOT(N2, V1) + d2;
    dv2 = _TRITRI_DOT(N2, V2) + d2;

#if _USE_EPSILON_TEST==TRUE
    if (Fabs(dv0) < _TRITRI_EPSILON) dv0 = ZERO;
    if (Fabs(dv1) < _TRITRI_EPSILON) dv1 = ZERO;
    if (Fabs(dv2) < _TRITRI_EPSILON) dv2 = ZERO;
#endif

    dv0dv1 = dv0 * dv1;
    dv0dv2 = dv0 * dv2;

    if (dv0dv1 > ZERO && dv0dv2 > ZERO) /* same sign on all of them + not equal 0 ? */
        return 0;                    /* no intersection occurs */

                                     /* compute direction of intersection line */
    _TRITRI_CROSS(D, N1, N2);

    /* compute and index to the largest component of D */
    max = Fabs(D[0]);
    index = 0;
    b = Fabs(D[1]);
    c = Fabs(D[2]);
    if (b > max) max = b, index = 1;
    if (c > max) max = c, index = 2;

    /* this is the simplified projection onto L*/
    vp0 = V0[index];
    vp1 = V1[index];
    vp2 = V2[index];

    up0 = U0[index];
    up1 = U1[index];
    up2 = U2[index];

    /* compute interval for triangle 1 */
    *coplanar = _compute_intervals_isectline(V0, V1, V2, vp0, vp1, vp2, dv0, dv1, dv2,
        dv0dv1, dv0dv2, &isect1[0], &isect1[1], isectpointA1, isectpointA2);
    if (*coplanar) return _coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);


    /* compute interval for triangle 2 */
    _compute_intervals_isectline(U0, U1, U2, up0, up1, up2, du0, du1, du2,
        du0du1, du0du2, &isect2[0], &isect2[1], isectpointB1, isectpointB2);

    _TRITRI_SORT2(isect1[0], isect1[1], smallest1);
    _TRITRI_SORT2(isect2[0], isect2[1], smallest2);

    if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;

    /* at this point, we know that the triangles intersect */

    if (isect2[0] < isect1[0])
    {
        if (smallest1 == 0) { _TRITRI_SET(isectpt1, isectpointA1); }
        else { _TRITRI_SET(isectpt1, isectpointA2); }

        if (isect2[1] < isect1[1])
        {
            if (smallest2 == 0) { _TRITRI_SET(isectpt2, isectpointB2); }
            else { _TRITRI_SET(isectpt2, isectpointB1); }
        }
        else
        {
            if (smallest1 == 0) { _TRITRI_SET(isectpt2, isectpointA2); }
            else { _TRITRI_SET(isectpt2, isectpointA1); }
        }
    }
    else
    {
        if (smallest2 == 0) { _TRITRI_SET(isectpt1, isectpointB1); }
        else { _TRITRI_SET(isectpt1, isectpointB2); }

        if (isect2[1] > isect1[1])
        {
            if (smallest1 == 0) { _TRITRI_SET(isectpt2, isectpointA2); }
            else { _TRITRI_SET(isectpt2, isectpointA1); }
        }
        else
        {
            if (smallest2 == 0) { _TRITRI_SET(isectpt2, isectpointB2); }
            else { _TRITRI_SET(isectpt2, isectpointB1); }
        }
    }
    return 1;
}
