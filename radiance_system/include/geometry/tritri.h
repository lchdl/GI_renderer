#pragma once
#include "basemath.h"
/*
reference:
https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/_TRITRI_isectline.txt
*/
/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 * updated: 2001-06-20 (added line of intersection)
 *
 * Here is a version withouts divisions (a little faster)
 * int NoDiv_TRITRIIsect(REAL V0[3],REAL V1[3],REAL V2[3],
 *                      REAL U0[3],REAL U1[3],REAL U2[3]);
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 * This version computes the line of intersection as well (if they are not coplanar):
 * int tri_tri_intersect_with_isectline(REAL V0[3],REAL V1[3],REAL V2[3],
 *				        REAL U0[3],REAL U1[3],REAL U2[3],int *coplanar,
 *				        REAL isectpt1[3],REAL isectpt2[3]);
 * coplanar returns whether the tris are coplanar
 * isectpt1, isectpt2 are the endpoints of the line of intersection
 */

 /*
Copyright 2020 Tomas Akenine-Möller

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and
to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <math.h>

/* if USE_EPSILON_TEST is true then we do a check:
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
#define _TRITRI_USE_EPSILON_TEST TRUE
#define _TRITRI_EPSILON REAL(1e-6)

/* some macros */
#define _TRITRI_CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 

#define _TRITRI_DOT(v1,v2)          (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define _TRITRI_SUB(dest,v1,v2)     dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2];
#define _TRITRI_ADD(dest,v1,v2)     dest[0]=v1[0]+v2[0]; dest[1]=v1[1]+v2[1]; dest[2]=v1[2]+v2[2];
#define _TRITRI_MULT(dest,v,factor) dest[0]=factor*v[0]; dest[1]=factor*v[1]; dest[2]=factor*v[2];
#define _TRITRI_SET(dest,src)       dest[0]=src[0]; dest[1]=src[1]; dest[2]=src[2];

/* sort so that a<=b */
#define _TRITRI_SORT(a,b) if(a>b) { REAL c; c=a; a=b; b=c; }

#define _TRITRI_ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);

#define _TRITRI_COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>ZERO)                                         \
  {                                                     \
    _ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>ZERO)                                    \
  {                                                     \
    _ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>ZERO || D0!=ZERO)                       \
  {                                                     \
    _ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=ZERO)                                     \
  {                                                     \
    _ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=ZERO)                                     \
  {                                                     \
    _ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    return _coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }

/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 
*/
#define _TRITRI_EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0) { if(e>=0 && e<=f) return 1;}             \
    else { if(e<=0 && e>=f) return 1; }}  

#define _TRITRI_EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  REAL Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  _TRITRI_EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  _TRITRI_EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  _TRITRI_EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define _TRITRI_POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  REAL a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

int _coplanar_tri_tri(REAL N[3], REAL V0[3], REAL V1[3], REAL V2[3],
    REAL U0[3], REAL U1[3], REAL U2[3]);

#define _TRITRI_NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
{ \
        if(D0D1>ZERO) \
        { \
                /* here we know that D0D2<=0.0 */ \
            /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
        } \
        else if(D0D2>ZERO)\
        { \
                /* here we know that d0d1<=0.0 */ \
            A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
        } \
        else if(D1*D2>ZERO || D0!=ZERO) \
        { \
                /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
                A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
        } \
        else if(D1!=ZERO) \
        { \
                A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
        } \
        else if(D2!=ZERO) \
        { \
                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
        } \
        else \
        { \
                /* triangles are coplanar */ \
                return _coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
        } \
}

int _tri_tri_intersect(REAL V0[3], REAL V1[3], REAL V2[3],
    REAL U0[3], REAL U1[3], REAL U2[3]);

/* sort so that a<=b */
#define _TRITRI_SORT2(a,b,smallest)       \
             if(a>b)       \
             {             \
               REAL c;    \
               c=a;        \
               a=b;        \
               b=c;        \
               smallest=1; \
             }             \
             else smallest=0;


void _isect2(REAL VTX0[3], REAL VTX1[3], REAL VTX2[3], REAL VV0, REAL VV1, REAL VV2,
    REAL D0, REAL D1, REAL D2, REAL *isect0, REAL *isect1, REAL isectpoint0[3], REAL isectpoint1[3]);

int _compute_intervals_isectline(REAL VERT0[3], REAL VERT1[3], REAL VERT2[3],
    REAL VV0, REAL VV1, REAL VV2, REAL D0, REAL D1, REAL D2,
    REAL D0D1, REAL D0D2, REAL *isect0, REAL *isect1,
    REAL isectpoint0[3], REAL isectpoint1[3]);

#define _TRITRI_COMPUTE_INTERVALS_ISECTLINE(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1,isectpoint0,isectpoint1) \
  if(D0D1>ZERO)                                         \
  {                                                     \
    _isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }

int _tri_tri_intersect_with_isectline(REAL V0[3], REAL V1[3], REAL V2[3],
    REAL U0[3], REAL U1[3], REAL U2[3], int *coplanar,
    REAL isectpt1[3], REAL isectpt2[3]);