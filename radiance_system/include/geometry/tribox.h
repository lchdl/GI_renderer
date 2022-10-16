#pragma once
#include "basemath.h"
/* 
Reference:
https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/
https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt 
*/
/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(REAL boxcenter[3],       */
/*          REAL boxhalfsize[3],REAL triverts[3][3]);   */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
/*
Copyright 2020 Tomas Akenine-Möller

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so, subject to the following 
conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
*/
#include <math.h>
#include <stdio.h>
#define _TRIBOX_X 0
#define _TRIBOX_Y 1
#define _TRIBOX_Z 2
#define _TRIBOX_CROSS(dest,v1,v2) \
dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define _TRIBOX_DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define _TRIBOX_SUB(dest,v1,v2) \
dest[0] = v1[0] - v2[0]; \
dest[1] = v1[1] - v2[1]; \
dest[2] = v1[2] - v2[2];

#define _TRIBOX_FINDMINMAX(x0,x1,x2,min,max) \
min = max = x0;                       \
if (x1 < min) min = x1;               \
    if (x1 > max) max = x1;           \
        if (x2 < min) min = x2;       \
            if (x2 > max) max = x2;

int _plane_box_overlap(REAL normal[3], REAL vert[3], REAL maxbox[3]);

/*======================== X-tests ========================*/
#define _TRIBOX_AXISTEST_X01(a, b, fa, fb)			                      \
p0 = a * v0[_TRIBOX_Y] - b * v0[_TRIBOX_Z];			       	                  \
p2 = a * v2[_TRIBOX_Y] - b * v2[_TRIBOX_Z];			       	                  \
if (p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; } \
rad = fa * boxhalfsize[_TRIBOX_Y] + fb * boxhalfsize[_TRIBOX_Z];                \
if (min > rad || max < -rad) return 0;

#define _TRIBOX_AXISTEST_X2(a, b, fa, fb)			                      \
p0 = a * v0[_TRIBOX_Y] - b * v0[_TRIBOX_Z];			                          \
p1 = a * v1[_TRIBOX_Y] - b * v1[_TRIBOX_Z];			       	                  \
if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; } \
rad = fa * boxhalfsize[_TRIBOX_Y] + fb * boxhalfsize[_TRIBOX_Z];                \
if (min > rad || max < -rad) return 0;

/*======================== Y-tests ========================*/
#define _TRIBOX_AXISTEST_Y02(a, b, fa, fb)			  	 	 	 	 	  \
p0 = -a * v0[_TRIBOX_X] + b * v0[_TRIBOX_Z];		      	  	 	 	 	 	  \
p2 = -a * v2[_TRIBOX_X] + b * v2[_TRIBOX_Z];	       	       	 	 	 	 	  \
if (p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; } \
rad = fa * boxhalfsize[_TRIBOX_X] + fb * boxhalfsize[_TRIBOX_Z];   	 	 	  \
if (min > rad || max < -rad) return 0;

#define _TRIBOX_AXISTEST_Y1(a, b, fa, fb)			   	 	 	 	 	  \
p0 = -a * v0[_TRIBOX_X] + b * v0[_TRIBOX_Z];		      		 	 	 	      \
p1 = -a * v1[_TRIBOX_X] + b * v1[_TRIBOX_Z];	     	       		 	 	 	  \
if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; } \
rad = fa * boxhalfsize[_TRIBOX_X] + fb * boxhalfsize[_TRIBOX_Z];   	 	 	  \
if (min > rad || max < -rad) return 0;

/*======================== Z-tests ========================*/
#define _TRIBOX_AXISTEST_Z12(a, b, fa, fb)			   	 	 	 	 	  \
p1 = a * v1[_TRIBOX_X] - b * v1[_TRIBOX_Y];			          	 	 	 	  \
p2 = a * v2[_TRIBOX_X] - b * v2[_TRIBOX_Y];			       	  	 	 	 	  \
if (p2 < p1) { min = p2; max = p1; } else { min = p1; max = p2; } \
rad = fa * boxhalfsize[_TRIBOX_X] + fb * boxhalfsize[_TRIBOX_Y];   	 	 	  \
if (min > rad || max < -rad) return 0;

#define _TRIBOX_AXISTEST_Z0(a, b, fa, fb)			   	 	 	 	 	  \
p0 = a * v0[_TRIBOX_X] - b * v0[_TRIBOX_Y];				  	 	 	 	 	  \
p1 = a * v1[_TRIBOX_X] - b * v1[_TRIBOX_Y];			          	 	 	 	  \
if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; } \
rad = fa * boxhalfsize[_TRIBOX_X] + fb * boxhalfsize[_TRIBOX_Y];   	 	      \
if (min > rad || max < -rad) return 0;

int _tri_box_overlap(REAL boxcenter[3], REAL boxhalfsize[3], REAL triverts[3][3]);
