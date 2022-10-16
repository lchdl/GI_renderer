#include "geometry/tribox.h"

int _plane_box_overlap(REAL normal[3], REAL vert[3], REAL maxbox[3])
{
    int q;
    REAL vmin[3], vmax[3], v;
    for (q = _TRIBOX_X; q <= _TRIBOX_Z; q++)
    {
        v = vert[q];
        if (normal[q] > ZERO)
        {
            vmin[q] = -maxbox[q] - v;
            vmax[q] = maxbox[q] - v;
        }
        else
        {
            vmin[q] = maxbox[q] - v;
            vmax[q] = -maxbox[q] - v;
        }
    }
    if (_TRIBOX_DOT(normal, vmin) > ZERO) return 0;
    if (_TRIBOX_DOT(normal, vmax) >= ZERO) return 1;
    return 0;
}

int _tri_box_overlap(REAL boxcenter[3], REAL boxhalfsize[3], REAL triverts[3][3])
{
    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */
    REAL v0[3], v1[3], v2[3];
    //   REAL axis[3];
    REAL min, max, p0, p1, p2, rad, fex, fey, fez;
    REAL normal[3], e0[3], e1[3], e2[3];
    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    _TRIBOX_SUB(v0, triverts[0], boxcenter);
    _TRIBOX_SUB(v1, triverts[1], boxcenter);
    _TRIBOX_SUB(v2, triverts[2], boxcenter);
    /* compute triangle edges */
    _TRIBOX_SUB(e0, v1, v0);      /* tri edge 0 */
    _TRIBOX_SUB(e1, v2, v1);      /* tri edge 1 */
    _TRIBOX_SUB(e2, v0, v2);      /* tri edge 2 */
    /* Bullet 3:  */
    /*  test the 9 tests first (this was faster) */
    fex = Fabs(e0[_TRIBOX_X]);
    fey = Fabs(e0[_TRIBOX_Y]);
    fez = Fabs(e0[_TRIBOX_Z]);
    _TRIBOX_AXISTEST_X01(e0[_TRIBOX_Z], e0[_TRIBOX_Y], fez, fey);
    _TRIBOX_AXISTEST_Y02(e0[_TRIBOX_Z], e0[_TRIBOX_X], fez, fex);
    _TRIBOX_AXISTEST_Z12(e0[_TRIBOX_Y], e0[_TRIBOX_X], fey, fex);
    fex = Fabs(e1[_TRIBOX_X]);
    fey = Fabs(e1[_TRIBOX_Y]);
    fez = Fabs(e1[_TRIBOX_Z]);
    _TRIBOX_AXISTEST_X01(e1[_TRIBOX_Z], e1[_TRIBOX_Y], fez, fey);
    _TRIBOX_AXISTEST_Y02(e1[_TRIBOX_Z], e1[_TRIBOX_X], fez, fex);
    _TRIBOX_AXISTEST_Z0(e1[_TRIBOX_Y], e1[_TRIBOX_X], fey, fex);
    fex = Fabs(e2[_TRIBOX_X]);
    fey = Fabs(e2[_TRIBOX_Y]);
    fez = Fabs(e2[_TRIBOX_Z]);
    _TRIBOX_AXISTEST_X2(e2[_TRIBOX_Z], e2[_TRIBOX_Y], fez, fey);
    _TRIBOX_AXISTEST_Y1(e2[_TRIBOX_Z], e2[_TRIBOX_X], fez, fex);
    _TRIBOX_AXISTEST_Z12(e2[_TRIBOX_Y], e2[_TRIBOX_X], fey, fex);
    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */
    /* test in X-direction */
    _TRIBOX_FINDMINMAX(v0[_TRIBOX_X], v1[_TRIBOX_X], v2[_TRIBOX_X], min, max);
    if (min > boxhalfsize[_TRIBOX_X] || max < -boxhalfsize[_TRIBOX_X]) return 0;
    /* test in Y-direction */
    _TRIBOX_FINDMINMAX(v0[_TRIBOX_Y], v1[_TRIBOX_Y], v2[_TRIBOX_Y], min, max);
    if (min > boxhalfsize[_TRIBOX_Y] || max < -boxhalfsize[_TRIBOX_Y]) return 0;
    /* test in Z-direction */
    _TRIBOX_FINDMINMAX(v0[_TRIBOX_Z], v1[_TRIBOX_Z], v2[_TRIBOX_Z], min, max);
    if (min > boxhalfsize[_TRIBOX_Z] || max < -boxhalfsize[_TRIBOX_Z]) return 0;
    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    _TRIBOX_CROSS(normal, e0, e1);
    if (!_plane_box_overlap(normal, v0, boxhalfsize)) return 0;
    return 1;   /* box and triangle overlaps */
}
