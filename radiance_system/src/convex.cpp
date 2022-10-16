#include "convex.h"

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/* structs for qsort */
struct _indexed_REAL { REAL val; int idx; };
struct _indexed_int { int val; int idx; };

/* internal functions definitions: */
int _convhull_cmp_asc_float(const void *a, const void *b) {
    _indexed_REAL *a1 = (_indexed_REAL*)a;
    _indexed_REAL *a2 = (_indexed_REAL*)b;
    if ((*a1).val < (*a2).val)return -1;
    else if ((*a1).val > (*a2).val)return 1;
    else return 0;
}
int _convhull_cmp_desc_float(const void *a, const void *b) {
    _indexed_REAL *a1 = (_indexed_REAL*)a;
    _indexed_REAL *a2 = (_indexed_REAL*)b;
    if ((*a1).val > (*a2).val)return -1;
    else if ((*a1).val < (*a2).val)return 1;
    else return 0;
}
int _convhull_cmp_asc_int(const void *a, const void *b) {
    _indexed_int *a1 = (_indexed_int*)a;
    _indexed_int *a2 = (_indexed_int*)b;
    if ((*a1).val < (*a2).val)return -1;
    else if ((*a1).val > (*a2).val)return 1;
    else return 0;
}
int _convhull_cmp_desc_int(const void *a, const void *b) {
    _indexed_int *a1 = (_indexed_int*)a;
    _indexed_int *a2 = (_indexed_int*)b;
    if ((*a1).val > (*a2).val)return -1;
    else if ((*a1).val < (*a2).val)return 1;
    else return 0;
}
void _convhull_sort_float(
    REAL* in_vec,      /* vector[len] to be sorted */
    REAL* out_vec,     /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices,   /* set to NULL if you don't need them */
    int len,           /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG    /* !1:ascending, 1:descending */
)
{
    int i;
    _indexed_REAL *data;

    data = (_indexed_REAL*)malloc(len * sizeof(_indexed_REAL));
    for (i = 0; i < len; i++) {
        data[i].val = in_vec[i];
        data[i].idx = i;
    }
    if (descendFLAG)
        qsort(data, len, sizeof(data[0]), _convhull_cmp_desc_float);
    else
        qsort(data, len, sizeof(data[0]), _convhull_cmp_asc_float);
    for (i = 0; i < len; i++) {
        if (out_vec != NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if (new_idices != NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}
void _convhull_sort_int
(
    int* in_vec,       /* vector[len] to be sorted */
    int* out_vec,      /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices,   /* set to NULL if you don't need them */
    int len,           /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG    /* !1:ascending, 1:descending */
)
{
    int i;
    _indexed_int *data;

    data = (_indexed_int*)malloc(len * sizeof(_indexed_int));
    for (i = 0; i < len; i++) {
        data[i].val = in_vec[i];
        data[i].idx = i;
    }
    if (descendFLAG)
        qsort(data, len, sizeof(data[0]), _convhull_cmp_desc_int);
    else
        qsort(data, len, sizeof(data[0]), _convhull_cmp_asc_int);
    for (i = 0; i < len; i++) {
        if (out_vec != NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if (new_idices != NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}
VEC3 _convhull_cross(VEC3* v1, VEC3* v2)
{
    VEC3 cross;
    cross.x = v1->y * v2->z - v1->z * v2->y;
    cross.y = v1->z * v2->x - v1->x * v2->z;
    cross.z = v1->x * v2->y - v1->y * v2->x;
    return cross;
}
REAL _convhull_det_4x4(REAL* m) {
    /* calculates the determinent of a 4x4 matrix */
    return
        m[3] * m[6] * m[9] * m[12] - m[2] * m[7] * m[9] * m[12] -
        m[3] * m[5] * m[10] * m[12] + m[1] * m[7] * m[10] * m[12] +
        m[2] * m[5] * m[11] * m[12] - m[1] * m[6] * m[11] * m[12] -
        m[3] * m[6] * m[8] * m[13] + m[2] * m[7] * m[8] * m[13] +
        m[3] * m[4] * m[10] * m[13] - m[0] * m[7] * m[10] * m[13] -
        m[2] * m[4] * m[11] * m[13] + m[0] * m[6] * m[11] * m[13] +
        m[3] * m[5] * m[8] * m[14] - m[1] * m[7] * m[8] * m[14] -
        m[3] * m[4] * m[9] * m[14] + m[0] * m[7] * m[9] * m[14] +
        m[1] * m[4] * m[11] * m[14] - m[0] * m[5] * m[11] * m[14] -
        m[2] * m[5] * m[8] * m[15] + m[1] * m[6] * m[8] * m[15] +
        m[2] * m[4] * m[9] * m[15] - m[0] * m[6] * m[9] * m[15] -
        m[1] * m[4] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];
}
void _convhull_create_sub_matrix
(/* Helper function for det_NxN()  */
    REAL* m,
    int N,
    int i,
    REAL* sub_m
)
{
    int j, k;
    for (j = N, k = 0; j < N * N; j++) {
        if (j % N != i) { /* i is the index to remove */
            sub_m[k] = m[j];
            k++;
        }
    }
}

REAL _convhull_det_NxN
(
    REAL* m,
    int d
)
{
    REAL sum;
    REAL sub_m[CONVHULL_ND_MAX_DIMENSIONS * CONVHULL_ND_MAX_DIMENSIONS];
    int sign;

    if (d == 0)
        return ONE;
    sum = ZERO;
    sign = 1;
    for (int i = 0; i < d; i++) {
        _convhull_create_sub_matrix(m, d, i, sub_m);
        sum += sign * m[i] * _convhull_det_NxN(sub_m, d - 1);
        sign *= -1;
    }
    return sum;
}

/* Calculates the coefficients of the equation of a PLANE in 3D.
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 */
void _convhull_plane_3d
(
    REAL* p,
    REAL* c,
    REAL* d
)
{
    int i, j, k, l;
    int r[3];
    REAL sign, det, norm_c;
    REAL pdiff[2][3], pdiff_s[2][2];

    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++)
            pdiff[i][j] = p[(i + 1) * 3 + j] - p[i * 3 + j];
    memset(c, 0, 3 * sizeof(REAL));
    sign = ONE;
    for (i = 0; i < 3; i++)
        r[i] = i;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 2; j++) {
            for (k = 0, l = 0; k < 3; k++) {
                if (r[k] != i) {
                    pdiff_s[j][l] = pdiff[j][k];
                    l++;
                }
            }
        }
        det = pdiff_s[0][0] * pdiff_s[1][1] - pdiff_s[1][0] * pdiff_s[0][1];
        c[i] = sign * det;
        sign *= -ONE;
    }
    norm_c = ZERO;
    for (i = 0; i < 3; i++)
        norm_c += (Pow(c[i], TWO));
    norm_c = Sqrt(norm_c);
    for (i = 0; i < 3; i++)
        c[i] /= norm_c;
    (*d) = ZERO;
    for (i = 0; i < 3; i++)
        (*d) += -p[i] * c[i];
}

/* Calculates the coefficients of the equation of a PLANE in ND.
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 */
void _convhull_plane_nd
(
    const int Nd,
    REAL* p,
    REAL* c,
    REAL* d
)
{
    int i, j, k, l;
    int r[CONVHULL_ND_MAX_DIMENSIONS];
    REAL sign, det, norm_c;
    REAL pdiff[CONVHULL_ND_MAX_DIMENSIONS - 1][CONVHULL_ND_MAX_DIMENSIONS], pdiff_s[(CONVHULL_ND_MAX_DIMENSIONS - 1)*(CONVHULL_ND_MAX_DIMENSIONS - 1)];

    if (Nd == 3) {
        _convhull_plane_3d(p, c, d);
        return;
    }

    for (i = 0; i < Nd - 1; i++)
        for (j = 0; j < Nd; j++)
            pdiff[i][j] = p[(i + 1)*Nd + j] - p[i*Nd + j];
    memset(c, 0, Nd * sizeof(REAL));
    sign = ONE;
    for (i = 0; i < Nd; i++)
        r[i] = i;
    for (i = 0; i < Nd; i++) {
        for (j = 0; j < Nd - 1; j++) {
            for (k = 0, l = 0; k < Nd; k++) {
                if (r[k] != i) {
                    pdiff_s[j*(Nd - 1) + l] = pdiff[j][k];
                    l++;
                }
            }
        }
        /* Determinant 1 dimension lower */
        if (Nd == 3)
            det = pdiff_s[0 * (Nd - 1) + 0] * pdiff_s[1 * (Nd - 1) + 1] - pdiff_s[1 * (Nd - 1) + 0] * pdiff_s[0 * (Nd - 1) + 1];
        else if (Nd == 5)
            det = _convhull_det_4x4((REAL*)pdiff_s);
        else {
            det = _convhull_det_NxN((REAL*)pdiff_s, Nd - 1);
        }
        c[i] = sign * det;
        sign *= -ONE;
    }
    norm_c = ZERO;
    for (i = 0; i < Nd; i++)
        norm_c += (Pow(c[i], TWO));
    norm_c = Sqrt(norm_c);
    for (i = 0; i < Nd; i++)
        c[i] /= norm_c;
    (*d) = ZERO;
    for (i = 0; i < Nd; i++)
        (*d) += -p[i] * c[i];
}

void _convhull_ismember
(
    int* pLeft,          /* left vector; nLeftElements x 1 */
    int* pRight,         /* right vector; nRightElements x 1 */
    int* pOut,           /* 0, unless pRight elements are present in pLeft then 1; nLeftElements x 1 */
    int nLeftElements,   /* number of elements in pLeft */
    int nRightElements   /* number of elements in pRight */
)
{
    int i, j;
    memset(pOut, 0, nLeftElements * sizeof(int));
    for (i = 0; i < nLeftElements; i++)
        for (j = 0; j < nRightElements; j++)
            if (pLeft[i] == pRight[j])
                pOut[i] = 1;
}

REAL _convhull_rnd(int x, int y)
{
    /*
    Reference(s):

    - Improvements to the canonical one-liner GLSL rand() for OpenGL ES 2.0
      http://byteblacksmith.com/improvements-to-the-canonical-one-liner-glsl-rand-for-opengl-es-2-0/
    */
    REAL a = REAL(12.9898);
    REAL b = REAL(78.233);
    REAL c = REAL(43758.5453);
    REAL dt = x * a + y * b;
    REAL sn = Fmod(dt, PI_LP);
    REAL intpart;
    return Modf(Sin(sn) * c, &intpart);
}



/* A C version of the 3D quickhull matlab implementation from here:
 * https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * (*out_faces) is returned as NULL, if triangulation fails *
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 * Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David P. Dobkin
 *             and Hannu Huhdanpaa, Geometry Center Technical Report GCG53, July 30, 1993"
 */

 /*
 builds the 3-D convexhull
 */
void _convhull_3d_build(
    /* input arguments */
    const VEC3* const in_vertices,      /* vector of input vertices; nVert x 1 */
    const int nVert,                    /* number of vertices */
    /* output arguments */
    int** out_faces,                    /* & of empty int*, output face indices; flat: nOut_faces x 3 */
    int* nOut_faces                     /* & of int, number of output face indices */
){
    int i, j, k, l, h;
    int nFaces, p, d;
    int* aVec, *faces;
    REAL dfi, v, max_p, min_p;
    REAL* points, *cf, *cfi, *df, *p_s, *span;

    if (nVert < 3 || in_vertices == NULL) {
        (*out_faces) = NULL;
        (*nOut_faces) = 0;
        return;
    }

    /* 3 dimensions. The code should theoretically work for >=2 dimensions, but "plane_3d" and "det_4x4" are hardcoded for 3,
     * so would need to be rewritten */
    d = 3;

    /* Add noise to the points */
    points = (REAL*)malloc(nVert*(d + 1) * sizeof(REAL));
    for (i = 0; i < nVert; i++) {
        for (j = 0; j < d; j++)
            points[i*(d + 1) + j] = in_vertices[i].e[j] + CONVHULL_NOISE_VALUE * _convhull_rnd(i, j); /* noise mitigates duplicates */
        points[i*(d + 1) + d] = ONE; /* add a last column of ones. Used only for determinant calculation */
    }
    /* Find the span */
    span = (REAL*)malloc(d * sizeof(REAL));
    for (j = 0; j < d; j++) {
        max_p = REAL(-2.23e+13); min_p = REAL(2.23e+13);
        for (i = 0; i < nVert; i++) {
            max_p = MAX(max_p, points[i*(d + 1) + j]);
            min_p = MIN(min_p, points[i*(d + 1) + j]);
        }
        span[j] = max_p - min_p;
        if (span[j] <= CONVHULL_NOISE_VALUE){
            /* If you hit this assertion error, then the input vertices do not span all 3 dimensions. Therefore the convex hull cannot be built.
             * In these cases, reduce the dimensionality of the points and call convhull_nd_build() instead with d<3 */
#ifdef CONVHULL_VERBOSE
            printf("ERROR: input vertices do not span all 3 dimensions. "
                   "Therefore the convex hull cannot be built. "
                   "In these cases, reduce the dimensionality of the "
                   "points and call convhull_nd_build() instead with d<3. "
                   "out_faces is set to NULL.\n");
#endif
            free(points);
            free(span);
            *out_faces = NULL;
            *nOut_faces = 0;
            return;
        }
    }

    /* The initial convex hull is a simplex with (d+1) facets, where d is the number of dimensions */
    nFaces = (d + 1);
    faces = (int*)calloc(nFaces*d, sizeof(int));
    aVec = (int*)malloc(nFaces * sizeof(int));
    for (i = 0; i < nFaces; i++)
        aVec[i] = i;

    /* Each column of cf contains the coefficients of a plane */
    cf = (REAL*)malloc(nFaces*d * sizeof(REAL));
    cfi = (REAL*)malloc(d * sizeof(REAL));
    df = (REAL*)malloc(nFaces * sizeof(REAL));
    p_s = (REAL*)malloc(d*d * sizeof(REAL));
    for (i = 0; i < nFaces; i++) {
        /* Set the indices of the points defining the face  */
        for (j = 0, k = 0; j < (d + 1); j++) {
            if (aVec[j] != i) {
                faces[i*d + k] = aVec[j];
                k++;
            }
        }

        /* Calculate and store the plane coefficients of the face */
        for (j = 0; j < d; j++)
            for (k = 0; k < d; k++)
                p_s[j*d + k] = points[(faces[i*d + j])*(d + 1) + k];

        /* Calculate and store the plane coefficients of the face */
        _convhull_plane_3d(p_s, cfi, &dfi);
        for (j = 0; j < d; j++)
            cf[i*d + j] = cfi[j];
        df[i] = dfi;
    }
    REAL *A;
    int *bVec, *fVec, *asfVec;
    int face_tmp[2];

    /* Check to make sure that faces are correctly oriented */
    bVec = (int*)malloc(4 * sizeof(int));
    for (i = 0; i < d + 1; i++)
        bVec[i] = i;

    /* A contains the coordinates of the points forming a simplex */
    A = (REAL*)calloc((d + 1)*(d + 1), sizeof(REAL));
    fVec = (int*)malloc((d + 1) * sizeof(int));
    asfVec = (int*)malloc((d + 1) * sizeof(int));
    for (k = 0; k < (d + 1); k++) {
        /* Get the point that is not on the current face (point p) */
        for (i = 0; i < d; i++)
            fVec[i] = faces[k*d + i];
        _convhull_sort_int(fVec, NULL, NULL, d, 0); /* sort accending */
        p = k;
        for (i = 0; i < d; i++)
            for (j = 0; j < (d + 1); j++)
                A[i*(d + 1) + j] = points[(faces[k*d + i])*(d + 1) + j];
        for (; i < (d + 1); i++)
            for (j = 0; j < (d + 1); j++)
                A[i*(d + 1) + j] = points[p*(d + 1) + j];

        /* det(A) determines the orientation of the face */
        v = _convhull_det_4x4(A);

        /* Orient so that each point on the original simplex can't see the opposite face */
        if (v < 0) {
            /* Reverse the order of the last two vertices to change the volume */
            for (j = 0; j < 2; j++)
                face_tmp[j] = faces[k*d + d - j - 1];
            for (j = 0; j < 2; j++)
                faces[k*d + d - j - 1] = face_tmp[1 - j];

            /* Modify the plane coefficients of the properly oriented faces */
            for (j = 0; j < d; j++)
                cf[k*d + j] = -cf[k*d + j];
            df[k] = -df[k];
            for (i = 0; i < d; i++)
                for (j = 0; j < (d + 1); j++)
                    A[i*(d + 1) + j] = points[(faces[k*d + i])*(d + 1) + j];
            for (; i < (d + 1); i++)
                for (j = 0; j < (d + 1); j++)
                    A[i*(d + 1) + j] = points[p*(d + 1) + j];
        }
    }

    /* Coordinates of the center of the point set */
    REAL* meanp, *absdist, *reldist, *desReldist;
    meanp = (REAL*)calloc(d, sizeof(REAL));
    for (i = d + 1; i < nVert; i++)
        for (j = 0; j < d; j++)
            meanp[j] += points[i*(d + 1) + j];
    for (j = 0; j < d; j++)
        meanp[j] = meanp[j] / (REAL)(nVert - d - 1);

    /* Absolute distance of points from the center */
    absdist = (REAL*)malloc((nVert - d - 1)*d * sizeof(REAL));
    for (i = d + 1, k = 0; i < nVert; i++, k++)
        for (j = 0; j < d; j++)
            absdist[k*d + j] = (points[i*(d + 1) + j] - meanp[j]) / span[j];

    /* Relative distance of points from the center */
    reldist = (REAL*)calloc((nVert - d - 1), sizeof(REAL));
    desReldist = (REAL*)malloc((nVert - d - 1) * sizeof(REAL));
    for (i = 0; i < (nVert - d - 1); i++)
        for (j = 0; j < d; j++)
            reldist[i] += Pow(absdist[i*d + j], TWO);

    /* Sort from maximum to minimum relative distance */
    int num_pleft, cnt;
    int* ind, *pleft;
    ind = (int*)malloc((nVert - d - 1) * sizeof(int));
    pleft = (int*)malloc((nVert - d - 1) * sizeof(int));
    _convhull_sort_float(reldist, desReldist, ind, (nVert - d - 1), 1);

    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first. */
    num_pleft = (nVert - d - 1);
    for (i = 0; i < num_pleft; i++)
        pleft[i] = ind[i] + d + 1;

    /* Loop over all remaining points that are not deleted. Deletion of points
     occurs every #iter2del# iterations of this while loop */
    memset(A, 0, (d + 1)*(d + 1) * sizeof(REAL));

    /* cnt is equal to the points having been selected without deletion of
     nonvisible points (i.e. points inside the current convex hull) */
    cnt = 0;

    /* The main loop for the quickhull algorithm */
    REAL detA;
    REAL* points_cf, *points_s;
    int* visible_ind, *visible, *nonvisible_faces, *f0, *face_s, *u, *gVec, *horizon, *hVec, *pp, *hVec_mem_face;
    int num_visible_ind, num_nonvisible_faces, n_newfaces, count, vis;
    int f0_sum, u_len, start, num_p, index, horizon_size1;
    int FUCKED;
    FUCKED = 0;
    u = horizon = NULL;
    nFaces = d + 1;
    visible_ind = (int*)malloc(nFaces * sizeof(int));
    points_cf = (REAL*)malloc(nFaces * sizeof(REAL));
    points_s = (REAL*)malloc(d * sizeof(REAL));
    face_s = (int*)malloc(d * sizeof(int));
    gVec = (int*)malloc(d * sizeof(int));
    while ((num_pleft > 0)) {
        /* i is the first point of the points left */
        i = pleft[0];

        /* Delete the point selected */
        for (j = 0; j < num_pleft - 1; j++)
            pleft[j] = pleft[j + 1];
        num_pleft--;
        if (num_pleft == 0)
            free(pleft);
        else
            pleft = (int*)realloc(pleft, num_pleft * sizeof(int));

        /* Update point selection counter */
        cnt++;

        /* find visible faces */
        for (j = 0; j < d; j++)
            points_s[j] = points[i*(d + 1) + j];
        points_cf = (REAL*)realloc(points_cf, nFaces * sizeof(REAL));
        visible_ind = (int*)realloc(visible_ind, nFaces * sizeof(int));
        for (j = 0; j < nFaces; j++) {
            points_cf[j] = 0;
            for (k = 0; k < d; k++)
                points_cf[j] += points_s[k] * cf[j*d + k];
        }
        num_visible_ind = 0;
        for (j = 0; j < nFaces; j++) {
            if (points_cf[j] + df[j] > ZERO) {
                num_visible_ind++; /* will sum to 0 if none are visible */
                visible_ind[j] = 1;
            }
            else
                visible_ind[j] = 0;
        }
        num_nonvisible_faces = nFaces - num_visible_ind;

        /* proceed if there are any visible faces */
        if (num_visible_ind != 0) {
            /* Find visible face indices */
            visible = (int*)malloc(num_visible_ind * sizeof(int));
            for (j = 0, k = 0; j < nFaces; j++) {
                if (visible_ind[j] == 1) {
                    visible[k] = j;
                    k++;
                }
            }

            /* Find nonvisible faces */
            nonvisible_faces = (int*)malloc(num_nonvisible_faces*d * sizeof(int));
            f0 = (int*)malloc(num_nonvisible_faces*d * sizeof(int));
            for (j = 0, k = 0; j < nFaces; j++) {
                if (visible_ind[j] == 0) {
                    for (l = 0; l < d; l++)
                        nonvisible_faces[k*d + l] = faces[j*d + l];
                    k++;
                }
            }

            /* Create horizon (count is the number of the edges of the horizon) */
            count = 0;
            for (j = 0; j < num_visible_ind; j++) {
                /* visible face */
                vis = visible[j];
                for (k = 0; k < d; k++)
                    face_s[k] = faces[vis*d + k];
                _convhull_sort_int(face_s, NULL, NULL, d, 0);
                _convhull_ismember(nonvisible_faces, face_s, f0, num_nonvisible_faces*d, d);
                u_len = 0;

                /* u are the nonvisible faces connected to the face v, if any */
                for (k = 0; k < num_nonvisible_faces; k++) {
                    f0_sum = 0;
                    for (l = 0; l < d; l++)
                        f0_sum += f0[k*d + l];
                    if (f0_sum == d - 1) {
                        u_len++;
                        if (u_len == 1)
                            u = (int*)malloc(u_len * sizeof(int));
                        else
                            u = (int*)realloc(u, u_len * sizeof(int));
                        u[u_len - 1] = k;
                    }
                }
                for (k = 0; k < u_len; k++) {
                    /* The boundary between the visible face v and the k(th) nonvisible face connected to the face v forms part of the horizon */
                    count++;
                    if (count == 1)
                        horizon = (int*)malloc(count*(d - 1) * sizeof(int));
                    else
                        horizon = (int*)realloc(horizon, count*(d - 1) * sizeof(int));
                    for (l = 0; l < d; l++)
                        gVec[l] = nonvisible_faces[u[k] * d + l];
                    for (l = 0, h = 0; l < d; l++) {
                        if (f0[u[k] * d + l]) {
                            horizon[(count - 1)*(d - 1) + h] = gVec[l];
                            h++;
                        }
                    }
                }
                if (u_len != 0)
                    free(u);
            }
            horizon_size1 = count;
            for (j = 0, l = 0; j < nFaces; j++) {
                if (!visible_ind[j]) {
                    /* Delete visible faces */
                    for (k = 0; k < d; k++)
                        faces[l*d + k] = faces[j*d + k];

                    /* Delete the corresponding plane coefficients of the faces */
                    for (k = 0; k < d; k++)
                        cf[l*d + k] = cf[j*d + k];
                    df[l] = df[j];
                    l++;
                }
            }

            /* Update the number of faces */
            nFaces = nFaces - num_visible_ind;
            faces = (int*)realloc(faces, nFaces*d * sizeof(int));
            cf = (REAL*)realloc(cf, nFaces*d * sizeof(REAL));
            df = (REAL*)realloc(df, nFaces * sizeof(REAL));

            /* start is the first row of the new faces */
            start = nFaces;

            /* Add faces connecting horizon to the new point */
            n_newfaces = horizon_size1;
            for (j = 0; j < n_newfaces; j++) {
                nFaces++;
                faces = (int*)realloc(faces, nFaces*d * sizeof(int));
                cf = (REAL*)realloc(cf, nFaces*d * sizeof(REAL));
                df = (REAL*)realloc(df, nFaces * sizeof(REAL));
                for (k = 0; k < d - 1; k++)
                    faces[(nFaces - 1)*d + k] = horizon[j*(d - 1) + k];
                faces[(nFaces - 1)*d + (d - 1)] = i;

                /* Calculate and store appropriately the plane coefficients of the faces */
                for (k = 0; k < d; k++)
                    for (l = 0; l < d; l++)
                        p_s[k*d + l] = points[(faces[(nFaces - 1)*d + k])*(d + 1) + l];
                _convhull_plane_3d(p_s, cfi, &dfi);
                for (k = 0; k < d; k++)
                    cf[(nFaces - 1)*d + k] = cfi[k];
                df[(nFaces - 1)] = dfi;
                if (nFaces > CONVHULL_MAX_NUM_FACES) {
#ifdef CONVHULL_VERBOSE
                    printf("ERROR: model is too complex, number of faces in convex hull "
                    "exceeds the maximum number of faces allowed. Please simplify the mesh "
                    "and try again!\n");
#endif
                    FUCKED = 1; /* hahaha */
                    nFaces = 0;
                    break;
                }
            }

            /* Orient each new face properly */
            hVec = (int*)malloc(nFaces * sizeof(int));
            hVec_mem_face = (int*)malloc(nFaces * sizeof(int));
            for (j = 0; j < nFaces; j++)
                hVec[j] = j;
            for (k = start; k < nFaces; k++) {
                for (j = 0; j < d; j++)
                    face_s[j] = faces[k*d + j];
                _convhull_sort_int(face_s, NULL, NULL, d, 0);
                _convhull_ismember(hVec, face_s, hVec_mem_face, nFaces, d);
                num_p = 0;
                for (j = 0; j < nFaces; j++)
                    if (!hVec_mem_face[j])
                        num_p++;
                pp = (int*)malloc(num_p * sizeof(int));
                for (j = 0, l = 0; j < nFaces; j++) {
                    if (!hVec_mem_face[j]) {
                        pp[l] = hVec[j];
                        l++;
                    }
                }
                index = 0;
                detA = ZERO;

                /* While new point is coplanar, choose another point */
                while (detA == ZERO) {
                    for (j = 0; j < d; j++)
                        for (l = 0; l < d + 1; l++)
                            A[j*(d + 1) + l] = points[(faces[k*d + j])*(d + 1) + l];
                    for (; j < d + 1; j++)
                        for (l = 0; l < d + 1; l++)
                            A[j*(d + 1) + l] = points[pp[index] * (d + 1) + l];
                    index++;
                    detA = _convhull_det_4x4(A);
                }

                /* Orient faces so that each point on the original simplex can't see the opposite face */
                if (detA < ZERO) {
                    /* If orientation is improper, reverse the order to change the volume sign */
                    for (j = 0; j < 2; j++)
                        face_tmp[j] = faces[k*d + d - j - 1];
                    for (j = 0; j < 2; j++)
                        faces[k*d + d - j - 1] = face_tmp[1 - j];

                    /* Modify the plane coefficients of the properly oriented faces */
                    for (j = 0; j < d; j++)
                        cf[k*d + j] = -cf[k*d + j];
                    df[k] = -df[k];
                    for (l = 0; l < d; l++)
                        for (j = 0; j < d + 1; j++)
                            A[l*(d + 1) + j] = points[(faces[k*d + l])*(d + 1) + j];
                    for (; l < d + 1; l++)
                        for (j = 0; j < d + 1; j++)
                            A[l*(d + 1) + j] = points[pp[index] * (d + 1) + j];
#ifdef CONVHULL_VERBOSE
                    /* Check */
                    detA = _convhull_det_4x4(A);
                    /* If you hit this assertion error, then the face cannot be properly orientated */
                    if (detA <= ZERO)
                        printf("WARNING: the face cannot be properly orientated.\n");
#endif
                }
                free(pp);
            }
            if (horizon_size1 > 0)
                free(horizon);
            free(f0);
            free(nonvisible_faces);
            free(visible);
            free(hVec);
            free(hVec_mem_face);
        }
        if (FUCKED) {
            break;
        }
    }

    /* output */
    if (FUCKED) {
        (*out_faces) = NULL;
        (*nOut_faces) = 0;
    }
    else {
        (*out_faces) = (int*)malloc(nFaces*d * sizeof(int));
        memcpy((*out_faces), faces, nFaces*d * sizeof(int));
        (*nOut_faces) = nFaces;
    }

    /* clean-up */
    free(visible_ind);
    free(points_cf);
    free(points_s);
    free(face_s);
    free(gVec);
    free(meanp);
    free(absdist);
    free(reldist);
    free(desReldist);
    free(ind);
    free(span);
    free(points);
    free(faces);
    free(aVec);
    free(cf);
    free(cfi);
    free(df);
    free(p_s);
    free(fVec);
    free(asfVec);
    free(bVec);
    free(A);
}

/*
exports the vertices, face indices, and face normals,
as an 'obj' file, ready for GPU (for 3d convexhulls only)
*/
void _convhull_3d_export_obj(
    /* input arguments */
    const VEC3* const vertices,         /* vector of input vertices; nVert x 1 */
    const int nVert,                    /* number of vertices */
    int* const faces,                   /* face indices; flat: nFaces x 3 */
    const int nFaces,                   /* number of faces in hull */
    const int keepOnlyUsedVerticesFLAG, /* 0: exports in_vertices, 1: exports only used vertices  */
    const char* obj_file                /* obj filename, WITH extension ".obj" */
){
    int i, j;
    FILE* fp;

    errno = 0;
    fp = fopen(obj_file, "w");

    if (fp == NULL) {
        printf("Error %d \n", errno);
        printf("It's null");
    }
    fprintf(fp, "o\n");
    REAL scale;
    VEC3 v1, v2, normal;

    /* export vertices */
    if (keepOnlyUsedVerticesFLAG) {
        for (i = 0; i < nFaces; i++)
            for (j = 0; j < 3; j++)
                fprintf(fp, "v %f %f %f\n", vertices[faces[i * 3 + j]].x,
                    vertices[faces[i * 3 + j]].y, vertices[faces[i * 3 + j]].z);
    }
    else {
        for (i = 0; i < nVert; i++)
            fprintf(fp, "v %f %f %f\n", vertices[i].x,
                vertices[i].y, vertices[i].z);
    }

    /* export the face normals */
    for (i = 0; i < nFaces; i++) {
        /* calculate cross product between v1-v0 and v2-v0 */
        v1 = vertices[faces[i * 3 + 1]];
        v2 = vertices[faces[i * 3 + 2]];
        v1.x -= vertices[faces[i * 3]].x;
        v1.y -= vertices[faces[i * 3]].y;
        v1.z -= vertices[faces[i * 3]].z;
        v2.x -= vertices[faces[i * 3]].x;
        v2.y -= vertices[faces[i * 3]].y;
        v2.z -= vertices[faces[i * 3]].z;
        normal = _convhull_cross(&v1, &v2);

        /* normalise to unit length */
        scale = ONE / (Sqrt(Pow(normal.x, TWO) + Pow(normal.y, TWO) + Pow(normal.z, TWO)) + (REAL)2.23e-9);
        normal.x *= scale;
        normal.y *= scale;
        normal.z *= scale;
        fprintf(fp, "vn %f %f %f\n", normal.x, normal.y, normal.z);
    }

    /* export the face indices */
    if (keepOnlyUsedVerticesFLAG) {
        for (i = 0; i < nFaces; i++) {
            /* vertices are in same order as the faces, and normals are in order */
            fprintf(fp, "f %u//%u %u//%u %u//%u\n",
                i * 3 + 1, i + 1,
                i * 3 + 1 + 1, i + 1,
                i * 3 + 2 + 1, i + 1);
        }
    }
    else {
        /* just normals are in order  */
        for (i = 0; i < nFaces; i++) {
            fprintf(fp, "f %u//%u %u//%u %u//%u\n",
                faces[i * 3] + 1, i + 1,
                faces[i * 3 + 1] + 1, i + 1,
                faces[i * 3 + 2] + 1, i + 1);
        }
    }
    fclose(fp);
}

/*
exports the vertices, face indices, and face normals, as an 'm' file,
for MatLab verification (for 3d convexhulls only)
*/
void _convhull_3d_export_m(
    /* input arguments */
    const VEC3* const vertices,         /* vector of input vertices; nVert x 1 */
    const int nVert,                    /* number of vertices */
    int* const faces,                   /* face indices; flat: nFaces x 3 */
    const int nFaces,                   /* number of faces in hull */
    char* const m_file                  /* m filename, WITH extension ".m" */
){
    int i;
    FILE* fp;
    fp = fopen(m_file, "w");

    /* save face indices and vertices for verification in matlab */
    fprintf(fp, "vertices = [\n");
    for (i = 0; i < nVert; i++)
        fprintf(fp, "%f, %f, %f;\n", vertices[i].x, vertices[i].y, vertices[i].z);
    fprintf(fp, "];\n\n\n");
    fprintf(fp, "faces = [\n");
    for (i = 0; i < nFaces; i++) {
        fprintf(fp, " %u, %u, %u;\n",
            faces[3 * i + 0] + 1,
            faces[3 * i + 1] + 1,
            faces[3 * i + 2] + 1);
    }
    fprintf(fp, "];\n\n\n");
    fclose(fp);
}

/*
reads an 'obj' file and extracts only the vertices
(for 3d convexhulls only)
*/
void _convhull_extract_vertices_from_obj(
    /* input arguments */
    const char* obj_file,               /* obj filename, WITH extension ".obj" */
    /* output arguments */
    VEC3** out_vertices,                /* & of empty VEC3*, output vertices; out_nVert x 1 */
    int* out_nVert                      /* & of int, number of vertices */
){
    FILE* fp;
    fp = fopen(obj_file, "r");

    /* determine number of vertices */
    unsigned int nVert = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        char* vexists = strstr(line, "v ");
        if (vexists != NULL)
            nVert++;
    }
    (*out_nVert) = nVert;
    (*out_vertices) = (VEC3*)malloc(nVert * sizeof(VEC3));

    /* extract the vertices */
    rewind(fp);
    int i = 0;
    int vertID, prev_char_isDigit, current_char_isDigit;
    char vert_char[256] = { 0 };
    while (fgets(line, sizeof(line), fp)) {
        char* vexists = strstr(line, "v ");
        if (vexists != NULL) {
            prev_char_isDigit = 0;
            vertID = -1;
            for (size_t j = 0; j < strlen(line) - 1; j++) {
                if (isdigit(line[j]) || line[j] == '.' || line[j] == '-' || line[j] == '+' || line[j] == 'E' || line[j] == 'e') {
                    vert_char[strlen(vert_char)] = line[j];
                    current_char_isDigit = 1;
                }
                else
                    current_char_isDigit = 0;
                if ((prev_char_isDigit && !current_char_isDigit) || j == strlen(line) - 2) {
                    vertID++;
                    if (vertID > 4) {
                        /* not a valid file */
                        free((*out_vertices));
                        (*out_vertices) = NULL;
                        (*out_nVert) = 0;
                        return;
                    }
                    (*out_vertices)[i].e[vertID] = (REAL)atof(vert_char);
                    memset(vert_char, 0, 256 * sizeof(char));
                }
                prev_char_isDigit = current_char_isDigit;
            }
            i++;
        }
    }
}


/* A C version of the ND quickhull matlab implementation from here:
 * https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * (*out_faces) is returned as NULL, if triangulation fails *
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 * Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David P. Dobkin
 *             and Hannu Huhdanpaa, Geometry Center Technical Report GCG53, July 30, 1993"
 */

/*
builds the N-Dimensional convexhull of a grid of points
*/
void _convhull_nd_build(
    /* input arguments */
    REAL* const in_points,              /* Matrix of points in 'd' dimensions; FLAT: nPoints x d */
    const int nPoints,                  /* number of points */
    const int d,                        /* Number of dimensions */
    /* output arguments */
    int** out_faces,                    /* (&) output face indices; FLAT: nOut_faces x d */
    REAL** out_cf,                      /* (&) contains the coefficients of the planes (set to NULL if not wanted);
                                           FLAT: nOut_faces x d */
    REAL** out_df,                      /* (&) contains the constant terms of the planes (set to NULL if not wanted);
                                           nOut_faces x 1 */
    int* nOut_faces                     /* (&) number of output face indices */
){
    int i, j, k, l, h;
    int nFaces, p;
    int* aVec, *faces;
    REAL dfi, v, max_p, min_p;
    REAL* points, *cf, *cfi, *df, *p_s, *span;

    /* Solution not possible... */
    if (nPoints <= d || in_points == NULL || d > CONVHULL_ND_MAX_DIMENSIONS) {
#ifdef CONVHULL_VERBOSE
        printf("ERROR: cannot build convex hull because a valid solution is not possible...\n"
               "Possible reasons are: 1) number of vertices is lower than dimension, 2) "
               "in_vertices is NULL, 3) dimension is too high (must lower than CONVHULL_ND_MAX_DIMENSIONS).\n");
#endif
        (*out_faces) = NULL;
        (*nOut_faces) = 0;
        if (out_cf != NULL) (*out_cf) = NULL;
        if (out_df != NULL) (*out_df) = NULL;
        return;
    }

    /* Add noise to the points */
    points = (REAL*)malloc(nPoints*(d + 1) * sizeof(REAL));
    for (i = 0; i < nPoints; i++) {
        for (j = 0; j < d; j++)
            points[i*(d + 1) + j] = in_points[i*d + j] + CONVHULL_NOISE_VALUE * _convhull_rnd(i, j);
        points[i*(d + 1) + d] = ONE; /* add a last column of ones. Used only for determinant calculation */
    }

    /* Find the span */
    span = (REAL*)malloc(d * sizeof(REAL));
    for (j = 0; j < d; j++) {
        max_p = REAL(-2.23e+13); min_p = REAL(2.23e+13);
        for (i = 0; i < nPoints; i++) {
            max_p = MAX(max_p, points[i*(d + 1) + j]);
            min_p = MIN(min_p, points[i*(d + 1) + j]);
        }
        span[j] = max_p - min_p;
        /* If you hit this assertion error, then the input vertices do not span all 'd' dimensions. Therefore the convex hull cannot be built.
         * In these cases, reduce the dimensionality of the points and call convhull_nd_build() instead with d<3 */
        if (span[j] <= CONVHULL_NOISE_VALUE) {
#ifdef CONVHULL_VERBOSE
            printf("ERROR: input vertices do not span all 'd' dimensions. "
                "Therefore the convex hull cannot be built. "
                "In these cases, reduce the dimensionality of the "
                "points and call convhull_nd_build() instead with d<3.\n");
#endif
            free(span);
            free(points);
            (*out_faces) = NULL;
            (*nOut_faces) = 0;
            if (out_cf != NULL) (*out_cf) = NULL;
            if (out_df != NULL) (*out_df) = NULL;
            return;
        }
    }

    /* The initial convex hull is a simplex with (d+1) facets, where d is the number of dimensions */
    nFaces = (d + 1);
    faces = (int*)calloc(nFaces*d, sizeof(int));
    aVec = (int*)malloc(nFaces * sizeof(int));
    for (i = 0; i < nFaces; i++)
        aVec[i] = i;

    /* Each column of cf contains the coefficients of a plane */
    cf = (REAL*)malloc(nFaces*d * sizeof(REAL));
    cfi = (REAL*)malloc(d * sizeof(REAL));
    df = (REAL*)malloc(nFaces * sizeof(REAL));
    p_s = (REAL*)malloc(d*d * sizeof(REAL));
    for (i = 0; i < nFaces; i++) {
        /* Set the indices of the points defining the face  */
        for (j = 0, k = 0; j < (d + 1); j++) {
            if (aVec[j] != i) {
                faces[i*d + k] = aVec[j];
                k++;
            }
        }

        /* Calculate and store the plane coefficients of the face */
        for (j = 0; j < d; j++)
            for (k = 0; k < d; k++)
                p_s[j*d + k] = points[(faces[i*d + j])*(d + 1) + k];

        /* Calculate and store the plane coefficients of the face */
        _convhull_plane_nd(d, p_s, cfi, &dfi);
        for (j = 0; j < d; j++)
            cf[i*d + j] = cfi[j];
        df[i] = dfi;
    }
    REAL *A;
    int *bVec, *fVec;
    int face_tmp[2];

    /* Check to make sure that faces are correctly oriented */
    bVec = (int*)malloc((d + 1) * sizeof(int));
    for (i = 0; i < d + 1; i++)
        bVec[i] = i;

    /* A contains the coordinates of the points forming a simplex */
    A = (REAL*)calloc((d + 1)*(d + 1), sizeof(REAL));
    fVec = (int*)malloc((d + 1) * sizeof(int));
    for (k = 0; k < (d + 1); k++) {
        /* Get the point that is not on the current face (point p) */
        for (i = 0; i < d; i++)
            fVec[i] = faces[k*d + i];
        _convhull_sort_int(fVec, NULL, NULL, d, 0); /* sort accending */
        p = k;
        for (i = 0; i < d; i++)
            for (j = 0; j < (d + 1); j++)
                A[i*(d + 1) + j] = points[(faces[k*d + i])*(d + 1) + j];
        for (; i < (d + 1); i++)
            for (j = 0; j < (d + 1); j++)
                A[i*(d + 1) + j] = points[p*(d + 1) + j];

        /* det(A) determines the orientation of the face */
        if (d == 3)
            v = _convhull_det_4x4(A);
        else
            v = _convhull_det_NxN(A, d + 1);

        /* Orient so that each point on the original simplex can't see the opposite face */
        if (v < 0) {
            /* Reverse the order of the last two vertices to change the volume */
            for (j = 0; j < 2; j++)
                face_tmp[j] = faces[k*d + d - j - 1];
            for (j = 0; j < 2; j++)
                faces[k*d + d - j - 1] = face_tmp[1 - j];

            /* Modify the plane coefficients of the properly oriented faces */
            for (j = 0; j < d; j++)
                cf[k*d + j] = -cf[k*d + j];
            df[k] = -df[k];
            for (i = 0; i < d; i++)
                for (j = 0; j < (d + 1); j++)
                    A[i*(d + 1) + j] = points[(faces[k*d + i])*(d + 1) + j];
            for (; i < (d + 1); i++)
                for (j = 0; j < (d + 1); j++)
                    A[i*(d + 1) + j] = points[p*(d + 1) + j];
        }
    }

    /* Coordinates of the center of the point set */
    REAL* meanp, *reldist, *desReldist, *absdist;
    meanp = (REAL*)calloc(d, sizeof(REAL));
    for (i = d + 1; i < nPoints; i++)
        for (j = 0; j < d; j++)
            meanp[j] += points[i*(d + 1) + j];
    for (j = 0; j < d; j++)
        meanp[j] = meanp[j] / (REAL)(nPoints - d - 1);

    /* Absolute distance of points from the center */
    absdist = (REAL*)malloc((nPoints - d - 1)*d * sizeof(REAL));
    for (i = d + 1, k = 0; i < nPoints; i++, k++)
        for (j = 0; j < d; j++)
            absdist[k*d + j] = (points[i*(d + 1) + j] - meanp[j]) / span[j];

    /* Relative distance of points from the center */
    reldist = (REAL*)calloc((nPoints - d - 1), sizeof(REAL));
    desReldist = (REAL*)malloc((nPoints - d - 1) * sizeof(REAL));
    for (i = 0; i < (nPoints - d - 1); i++)
        for (j = 0; j < d; j++)
            reldist[i] += Pow(absdist[i*d + j], TWO);

    /* Sort from maximum to minimum relative distance */
    int num_pleft, cnt;
    int* ind, *pleft;
    ind = (int*)malloc((nPoints - d - 1) * sizeof(int));
    pleft = (int*)malloc((nPoints - d - 1) * sizeof(int));
    _convhull_sort_float(reldist, desReldist, ind, (nPoints - d - 1), 1);

    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first. */
    num_pleft = (nPoints - d - 1);
    for (i = 0; i < num_pleft; i++)
        pleft[i] = ind[i] + d + 1;

    /* Loop over all remaining points that are not deleted. Deletion of points
     occurs every #iter2del# iterations of this while loop */
    memset(A, 0, (d + 1)*(d + 1) * sizeof(REAL));

    /* cnt is equal to the points having been selected without deletion of
     nonvisible points (i.e. points inside the current convex hull) */
    cnt = 0;

    /* The main loop for the quickhull algorithm */
    REAL detA;
    REAL* points_cf, *points_s;
    int* visible_ind, *visible, *nonvisible_faces, *f0, *face_s, *u, *gVec, *horizon, *hVec, *pp, *hVec_mem_face;
    int num_visible_ind, num_nonvisible_faces, n_newfaces, count, vis;
    int f0_sum, u_len, start, num_p, index, horizon_size1;
    int FUCKED;
    FUCKED = 0;
    u = horizon = NULL;
    nFaces = d + 1;
    visible_ind = (int*)malloc(nFaces * sizeof(int));
    points_cf = (REAL*)malloc(nFaces * sizeof(REAL));
    points_s = (REAL*)malloc(d * sizeof(REAL));
    face_s = (int*)malloc(d * sizeof(int));
    gVec = (int*)malloc(d * sizeof(int));
    while ((num_pleft > 0)) {
        /* i is the first point of the points left */
        i = pleft[0];

        /* Delete the point selected */
        for (j = 0; j < num_pleft - 1; j++)
            pleft[j] = pleft[j + 1];
        num_pleft--;
        if (num_pleft == 0)
            free(pleft);
        else
            pleft = (int*)realloc(pleft, num_pleft * sizeof(int));

        /* Update point selection counter */
        cnt++;

        /* find visible faces */
        for (j = 0; j < d; j++)
            points_s[j] = points[i*(d + 1) + j];
        points_cf = (REAL*)realloc(points_cf, nFaces * sizeof(REAL));
        visible_ind = (int*)realloc(visible_ind, nFaces * sizeof(int));
        for (j = 0; j < nFaces; j++) {
            points_cf[j] = 0;
            for (k = 0; k < d; k++)
                points_cf[j] += points_s[k] * cf[j*d + k];
        }
        num_visible_ind = 0;
        for (j = 0; j < nFaces; j++) {
            if (points_cf[j] + df[j] > ZERO) {
                num_visible_ind++; /* will sum to 0 if none are visible */
                visible_ind[j] = 1;
            }
            else
                visible_ind[j] = 0;
        }
        num_nonvisible_faces = nFaces - num_visible_ind;

        /* proceed if there are any visible faces */
        if (num_visible_ind != 0) {
            /* Find visible face indices */
            visible = (int*)malloc(num_visible_ind * sizeof(int));
            for (j = 0, k = 0; j < nFaces; j++) {
                if (visible_ind[j] == 1) {
                    visible[k] = j;
                    k++;
                }
            }

            /* Find nonvisible faces */
            nonvisible_faces = (int*)malloc(num_nonvisible_faces*d * sizeof(int));
            f0 = (int*)malloc(num_nonvisible_faces*d * sizeof(int));
            for (j = 0, k = 0; j < nFaces; j++) {
                if (visible_ind[j] == 0) {
                    for (l = 0; l < d; l++)
                        nonvisible_faces[k*d + l] = faces[j*d + l];
                    k++;
                }
            }

            /* Create horizon (count is the number of the edges of the horizon) */
            count = 0;
            for (j = 0; j < num_visible_ind; j++) {
                /* visible face */
                vis = visible[j];
                for (k = 0; k < d; k++)
                    face_s[k] = faces[vis*d + k];
                _convhull_sort_int(face_s, NULL, NULL, d, 0);
                _convhull_ismember(nonvisible_faces, face_s, f0, num_nonvisible_faces*d, d);
                u_len = 0;

                /* u are the nonvisible faces connected to the face v, if any */
                for (k = 0; k < num_nonvisible_faces; k++) {
                    f0_sum = 0;
                    for (l = 0; l < d; l++)
                        f0_sum += f0[k*d + l];
                    if (f0_sum == d - 1) {
                        u_len++;
                        if (u_len == 1)
                            u = (int*)malloc(u_len * sizeof(int));
                        else
                            u = (int*)realloc(u, u_len * sizeof(int));
                        u[u_len - 1] = k;
                    }
                }
                for (k = 0; k < u_len; k++) {
                    /* The boundary between the visible face v and the k(th) nonvisible face connected to the face v forms part of the horizon */
                    count++;
                    if (count == 1)
                        horizon = (int*)malloc(count*(d - 1) * sizeof(int));
                    else
                        horizon = (int*)realloc(horizon, count*(d - 1) * sizeof(int));
                    for (l = 0; l < d; l++)
                        gVec[l] = nonvisible_faces[u[k] * d + l];
                    for (l = 0, h = 0; l < d; l++) {
                        if (f0[u[k] * d + l]) {
                            horizon[(count - 1)*(d - 1) + h] = gVec[l];
                            h++;
                        }
                    }
                }
                if (u_len != 0)
                    free(u);
            }
            horizon_size1 = count;
            for (j = 0, l = 0; j < nFaces; j++) {
                if (!visible_ind[j]) {
                    /* Delete visible faces */
                    for (k = 0; k < d; k++)
                        faces[l*d + k] = faces[j*d + k];

                    /* Delete the corresponding plane coefficients of the faces */
                    for (k = 0; k < d; k++)
                        cf[l*d + k] = cf[j*d + k];
                    df[l] = df[j];
                    l++;
                }
            }

            /* Update the number of faces */
            nFaces = nFaces - num_visible_ind;
            faces = (int*)realloc(faces, nFaces*d * sizeof(int));
            cf = (REAL*)realloc(cf, nFaces*d * sizeof(REAL));
            df = (REAL*)realloc(df, nFaces * sizeof(REAL));

            /* start is the first row of the new faces */
            start = nFaces;

            /* Add faces connecting horizon to the new point */
            n_newfaces = horizon_size1;
            for (j = 0; j < n_newfaces; j++) {
                nFaces++;
                faces = (int*)realloc(faces, nFaces*d * sizeof(int));
                cf = (REAL*)realloc(cf, nFaces*d * sizeof(REAL));
                df = (REAL*)realloc(df, nFaces * sizeof(REAL));
                for (k = 0; k < d - 1; k++)
                    faces[(nFaces - 1)*d + k] = horizon[j*(d - 1) + k];
                faces[(nFaces - 1)*d + (d - 1)] = i;

                /* Calculate and store appropriately the plane coefficients of the faces */
                for (k = 0; k < d; k++)
                    for (l = 0; l < d; l++)
                        p_s[k*d + l] = points[(faces[(nFaces - 1)*d + k])*(d + 1) + l];
                _convhull_plane_nd(d, p_s, cfi, &dfi);
                for (k = 0; k < d; k++)
                    cf[(nFaces - 1)*d + k] = cfi[k];
                df[(nFaces - 1)] = dfi;
                if (nFaces > CONVHULL_MAX_NUM_FACES) {
                    FUCKED = 1;
                    nFaces = 0;
                    break;
                }
            }

            /* Orient each new face properly */
            hVec = (int*)malloc(nFaces * sizeof(int));
            hVec_mem_face = (int*)malloc(nFaces * sizeof(int));
            for (j = 0; j < nFaces; j++)
                hVec[j] = j;
            for (k = start; k < nFaces; k++) {
                for (j = 0; j < d; j++)
                    face_s[j] = faces[k*d + j];
                _convhull_sort_int(face_s, NULL, NULL, d, 0);
                _convhull_ismember(hVec, face_s, hVec_mem_face, nFaces, d);
                num_p = 0;
                for (j = 0; j < nFaces; j++)
                    if (!hVec_mem_face[j])
                        num_p++;
                pp = (int*)malloc(num_p * sizeof(int));
                for (j = 0, l = 0; j < nFaces; j++) {
                    if (!hVec_mem_face[j]) {
                        pp[l] = hVec[j];
                        l++;
                    }
                }
                index = 0;
                detA = ZERO;

                /* While new point is coplanar, choose another point */
                while (detA == ZERO) {
                    for (j = 0; j < d; j++)
                        for (l = 0; l < d + 1; l++)
                            A[j*(d + 1) + l] = points[(faces[k*d + j])*(d + 1) + l];
                    for (; j < d + 1; j++)
                        for (l = 0; l < d + 1; l++)
                            A[j*(d + 1) + l] = points[pp[index] * (d + 1) + l];
                    index++;
                    if (d == 3)
                        detA = _convhull_det_4x4(A);
                    else
                        detA = _convhull_det_NxN((REAL*)A, d + 1);
                }

                /* Orient faces so that each point on the original simplex can't see the opposite face */
                if (detA < ZERO) {
                    /* If orientation is improper, reverse the order to change the volume sign */
                    for (j = 0; j < 2; j++)
                        face_tmp[j] = faces[k*d + d - j - 1];
                    for (j = 0; j < 2; j++)
                        faces[k*d + d - j - 1] = face_tmp[1 - j];

                    /* Modify the plane coefficients of the properly oriented faces */
                    for (j = 0; j < d; j++)
                        cf[k*d + j] = -cf[k*d + j];
                    df[k] = -df[k];
                    for (l = 0; l < d; l++)
                        for (j = 0; j < d + 1; j++)
                            A[l*(d + 1) + j] = points[(faces[k*d + l])*(d + 1) + j];
                    for (; l < d + 1; l++)
                        for (j = 0; j < d + 1; j++)
                            A[l*(d + 1) + j] = points[pp[index] * (d + 1) + j];
#ifdef CONVHULL_VERBOSE
                    /* Check */
                    if (d == 3)
                        detA = _convhull_det_4x4(A);
                    else
                        detA = _convhull_det_NxN((REAL*)A, d + 1);
                    /* If you hit this assertion error, then the face cannot be properly orientated and building the convex hull is likely impossible */
                    if (detA <= ZERO)
                        printf("WARNING: the face cannot be properly orientated and building the convex hull is likely impossible.\n");
#endif
                }
                free(pp);
            }
            if (horizon_size1 > 0)
                free(horizon);
            free(f0);
            free(nonvisible_faces);
            free(visible);
            free(hVec);
            free(hVec_mem_face);
        }
        if (FUCKED) {
            break;
        }
    }

    /* output */
    if (FUCKED) {
        (*out_faces) = NULL;
        if (out_cf != NULL) (*out_cf) = NULL;
        if (out_df != NULL) (*out_df) = NULL;
        (*nOut_faces) = 0;
    }
    else {
        (*out_faces) = (int*)malloc(nFaces*d * sizeof(int));
        memcpy((*out_faces), faces, nFaces*d * sizeof(int));
        (*nOut_faces) = nFaces;
        if (out_cf != NULL) {
            (*out_cf) = (REAL*)malloc(nFaces*d * sizeof(REAL));
            memcpy((*out_cf), cf, nFaces*d * sizeof(REAL));
        }
        if (out_df != NULL) {
            (*out_df) = (REAL*)malloc(nFaces * sizeof(REAL));
            memcpy((*out_df), df, nFaces * sizeof(REAL));
        }
    }

    /* clean-up */
    free(visible_ind);
    free(points_cf);
    free(points_s);
    free(face_s);
    free(gVec);
    free(meanp);
    free(absdist);
    free(reldist);
    free(desReldist);
    free(ind);
    free(span);
    free(points);
    free(faces);
    free(aVec);
    free(cf);
    free(cfi);
    free(df);
    free(p_s);
    free(fVec);
    free(bVec);
    free(A);
}


/*
Computes the Delaunay triangulation (mesh) of an arrangement
of points in N-dimensional space
*/
void _convhull_delaunay_nd_mesh(
    /* input Arguments */
    const REAL* points,                 /* The input points; FLAT: nPoints x nd */
    const int nPoints,                  /* Number of points */
    const int nd,                       /* The number of dimensions */
    /* output Arguments */
    int** Mesh,                         /* (&) the indices defining the Delaunay triangulation of the points;
                                           FLAT: nMesh x (nd+1) */
    int* nMesh                          /* (&) Number of triangulations */
){
    int i, j, k, nHullFaces, maxW_idx, nVisible;
    int* hullfaces;
    REAL w0, w_optimal, w_optimal2;
    REAL* projpoints, *cf, *df, *p0, *p, *visible;

    /* Project the N-dimensional points onto a N+1-dimensional paraboloid */
    projpoints = (REAL*)malloc(nPoints*(nd + 1) * sizeof(REAL));
    for (i = 0; i < nPoints; i++) {
        projpoints[i*(nd + 1) + nd] = ZERO;
        for (j = 0; j < nd; j++) {
            projpoints[i*(nd + 1) + j] = (REAL)points[i*nd + j] + CONVHULL_NOISE_VALUE * _convhull_rnd(i, j);
            projpoints[i*(nd + 1) + nd] += (projpoints[i*(nd + 1) + j] * projpoints[i*(nd + 1) + j]); /* w vector */
        }
    }

    /* The N-dimensional delaunay triangulation requires first computing the convex hull of this N+1-dimensional paraboloid */
    hullfaces = NULL;
    cf = df = NULL;
    _convhull_nd_build(projpoints, nPoints, nd + 1, &hullfaces, &cf, &df, &nHullFaces);

    /* Find the coordinates of the point with the maximum (N+1 dimension) coordinate (i.e. the w vector) */
    REAL maxVal;
    maxVal = REAL(-2.23e13);
    maxW_idx = -1;
    for (i = 0; i < nPoints; i++) {
        if (projpoints[i*(nd + 1) + nd] > maxVal) {
            maxVal = projpoints[i*(nd + 1) + nd];
            maxW_idx = i;
        }
    }
    if (maxW_idx == -1) {
#ifdef CONVHULL_VERBOSE
        printf("ERROR: maxW_idx==-1.\n");
#endif
        free(projpoints);
        *Mesh = NULL;
        *nMesh = 0;
        return;
    }
    w0 = projpoints[maxW_idx*(nd + 1) + nd];
    p0 = (REAL*)malloc(nd * sizeof(REAL));
    for (j = 0; j < nd; j++)
        p0[j] = projpoints[maxW_idx*(nd + 1) + j];

    /* Find the point where the plane tangent to the point (p0,w0) on the paraboloid crosses the w axis.
     * This is the point that can see the entire lower hull. */
    w_optimal = ZERO;
    for (j = 0; j < nd; j++)
        w_optimal += (TWO*Pow(p0[j], TWO));
    w_optimal = w0 - w_optimal;

    /* Subtract 1000 times the absolute value of w_optimal to ensure that the point where the tangent plane
     * crosses the w axis will see all points on the lower hull. This avoids numerical roundoff errors. */
    w_optimal2 = w_optimal - REAL(1000.0) * Fabs(w_optimal);

    /* Set the point where the tangent plane crosses the w axis */
    p = (REAL*)calloc((nd + 1), sizeof(REAL));
    p[nd] = w_optimal2;

    /* Find all faces that are visible from this point */
    visible = (REAL*)malloc(nHullFaces * sizeof(REAL));
    for (i = 0; i < nHullFaces; i++) {
        visible[i] = ZERO;
        for (j = 0; j < nd + 1; j++)
            visible[i] += cf[i*(nd + 1) + j] * p[j];
    }
    nVisible = 0;
    for (j = 0; j < nHullFaces; j++) {
        visible[j] += df[j];
        if (visible[j] > ZERO)
            nVisible++;
    }

    /* Output */
    (*nMesh) = nVisible;
    if (nVisible > 0) {
        (*Mesh) = (int*)malloc(nVisible*(nd + 1) * sizeof(int));
        for (i = 0, j = 0; i < nHullFaces; i++) {
            if (visible[i] > ZERO) {
                for (k = 0; k < nd + 1; k++)
                    (*Mesh)[j*(nd + 1) + k] = hullfaces[i*(nd + 1) + k];
                j++;
            }
        }
#ifdef CONVHULL_VERBOSE
        if (j != nVisible) {
            printf("WAARNING: j!=nVisible .\n");
        }
#endif
    }

    /* clean up */
    free(projpoints);
    free(hullfaces);
    free(cf);
    free(df);
    free(p0);
    free(p);
    free(visible);
}

convex buildConvex3d(Array<VEC3>& vertices)
{
    convex cvex;
    int* faceIdxs = NULL;
    int nFaces = 0;
    _convhull_3d_build(vertices.data(), vertices.size(), &faceIdxs, &nFaces);
    if (faceIdxs == NULL) /* cannot build convex hull, return empty shape */
        return convex();

    /* pack all vertices */
    Array<int> vMap; /* packed vertex index -> unpacked index */
    for (int iFace = 0; iFace < nFaces; iFace++) {
        for (int iTriangle = 0; iTriangle < 3; iTriangle++) {
            int unpackedIndex, inArray = 0;
            unpackedIndex = faceIdxs[iFace * 3 + iTriangle];
            for (int j = 0; j < vMap.size(); j++) {
                if (vMap[j] == unpackedIndex) { /* this vertex is already in array */
                    inArray = 1; break;
                }
            }
            if (!inArray) vMap.append(unpackedIndex);
        }
    }
    /* fill in pHull */
    for (int i = 0; i < vMap.size(); i++) cvex.vertices.append(vertices[vMap[i]]);
    /* pack face indices */
    int* faceIdxsPacked = (int*)malloc(nFaces * 3 * sizeof(int));
    for (int i = 0; i < vMap.size(); i++) {
        for (int j = 0; j < nFaces * 3; j++) {
            if (faceIdxs[j] == vMap[i]) faceIdxsPacked[j] = i;
        }
    }
    /* now vMap is useless, clear it */
    vMap.clear();
    /* fill_ in face indices */
    for (int iF = 0; iF < nFaces; iF++) {
        INT3 f;
        f.x = faceIdxsPacked[iF * 3 + 0];
        f.y = faceIdxsPacked[iF * 3 + 1];
        f.z = faceIdxsPacked[iF * 3 + 2];
        cvex.faces.append(f);
    }

    free(faceIdxs);
    free(faceIdxsPacked);

    return cvex;
}
