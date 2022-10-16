/*

convex.h : extract convex hull shape from a point cloud.
Original source code is from:

https://github.com/leomccormack/convhull_3d/blob/master/convhull_3d.h

*/

#pragma once
#include "basemath.h"
#include "basedefs.h"

#define CONVHULL_MAX_NUM_FACES      100000
#define CONVHULL_ND_MAX_DIMENSIONS  5
#ifdef USE_DOUBLE_PRECISION
#define CONVHULL_NOISE_VALUE        double(0.0000001)
#else
#define CONVHULL_NOISE_VALUE        float(0.00001)
#endif
#define CONVHULL_VERBOSE /* enable some verbose outputs when failure */

/* definition of a 3D convex shape */
struct convex {
    Array<VEC3> vertices; /* all vertices of the convex hull shape */
    Array<INT3> faces; /* face indices, each INT3 represents all vertex indices of a triangle */
};

convex buildConvex3d(Array<VEC3>& vertices); /* build convex shape from point cloud */

