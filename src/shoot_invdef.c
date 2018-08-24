/*
 * $Id: shoot_invdef.c 7408 2018-08-24 14:54:57Z john $
 * John Ashburner
 */

/*  Invert a deformation field.

    The field is assumed to consist of a piecewise affine transformations, whereby
    each cube jointing 8 neighbouring voxels contains 8 tetrahedra.  The mapping
    within each tetrahedron is assumed to be affine.
    More documentation can be found in the appendix of:
    J. Ashburner, J. Andersson and K. J. Friston (2000).
    "Image Registration using a Symmetric Prior - in Three-Dimensions".
    Human Brain Mapping 9(4):212-225.
*/

#include <math.h>
#include "mex.h"

#define MAXV 16384

static void invertX(float X[][3], float IX[][4])
/* X is a matrix containing the co-ordinates of the four vertices of a tetrahedron.
   IX = inv([X ; 1 1 1 1]);  */
{
    float id;
    id = X[0][0]*(X[3][1]*(X[1][2]-X[2][2])+X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2]))+
         X[1][0]*(X[3][2]*(X[0][1]-X[2][1])+X[0][2]*(X[2][1]-X[3][1])+X[2][2]*(X[3][1]-X[0][1]))+
         X[2][0]*(X[0][1]*(X[1][2]-X[3][2])+X[3][1]*(X[0][2]-X[1][2])+X[1][1]*(X[3][2]-X[0][2]))+
         X[3][0]*(X[1][2]*(X[2][1]-X[0][1])+X[0][2]*(X[1][1]-X[2][1])+X[2][2]*(X[0][1]-X[1][1]));
    id = 1.0f/id;
    IX[0][0] = id*(X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2])+X[3][1]*(X[1][2]-X[2][2]));
    IX[0][1] = id*(X[0][1]*(X[3][2]-X[2][2])+X[2][1]*(X[0][2]-X[3][2])+X[3][1]*(X[2][2]-X[0][2]));
    IX[0][2] = id*(X[0][1]*(X[1][2]-X[3][2])+X[1][1]*(X[3][2]-X[0][2])+X[3][1]*(X[0][2]-X[1][2]));
    IX[0][3] = id*(X[0][1]*(X[2][2]-X[1][2])+X[1][1]*(X[0][2]-X[2][2])+X[2][1]*(X[1][2]-X[0][2]));
    IX[1][0] = id*(X[1][0]*(X[3][2]-X[2][2])+X[2][0]*(X[1][2]-X[3][2])+X[3][0]*(X[2][2]-X[1][2]));
    IX[1][1] = id*(X[0][0]*(X[2][2]-X[3][2])+X[2][0]*(X[3][2]-X[0][2])+X[3][0]*(X[0][2]-X[2][2]));
    IX[1][2] = id*(X[0][0]*(X[3][2]-X[1][2])+X[1][0]*(X[0][2]-X[3][2])+X[3][0]*(X[1][2]-X[0][2]));
    IX[1][3] = id*(X[0][0]*(X[1][2]-X[2][2])+X[1][0]*(X[2][2]-X[0][2])+X[2][0]*(X[0][2]-X[1][2]));
    IX[2][0] = id*(X[1][0]*(X[2][1]-X[3][1])+X[2][0]*(X[3][1]-X[1][1])+X[3][0]*(X[1][1]-X[2][1]));
    IX[2][1] = id*(X[0][0]*(X[3][1]-X[2][1])+X[2][0]*(X[0][1]-X[3][1])+X[3][0]*(X[2][1]-X[0][1]));
    IX[2][2] = id*(X[0][0]*(X[1][1]-X[3][1])+X[1][0]*(X[3][1]-X[0][1])+X[3][0]*(X[0][1]-X[1][1]));
    IX[2][3] = id*(X[0][0]*(X[2][1]-X[1][1])+X[1][0]*(X[0][1]-X[2][1])+X[2][0]*(X[1][1]-X[0][1]));
    IX[3][0] = id*(X[1][0]*(X[2][2]*X[3][1]-X[3][2]*X[2][1])+
                   X[2][0]*(X[3][2]*X[1][1]-X[1][2]*X[3][1])+
               X[3][0]*(X[1][2]*X[2][1]-X[2][2]*X[1][1]));
    IX[3][1] = id*(X[0][0]*(X[3][2]*X[2][1]-X[2][2]*X[3][1])+
               X[2][0]*(X[0][2]*X[3][1]-X[3][2]*X[0][1])+
               X[3][0]*(X[2][2]*X[0][1]-X[0][2]*X[2][1]));
    IX[3][2] = id*(X[0][0]*(X[1][2]*X[3][1]-X[3][2]*X[1][1])+
               X[1][0]*(X[3][2]*X[0][1]-X[0][2]*X[3][1])+
               X[3][0]*(X[0][2]*X[1][1]-X[1][2]*X[0][1]));
    IX[3][3] = id*(X[0][0]*(X[2][2]*X[1][1]-X[1][2]*X[2][1])+
               X[1][0]*(X[0][2]*X[2][1]-X[2][2]*X[0][1])+
               X[2][0]*(X[1][2]*X[0][1]-X[0][2]*X[1][1]));
}

static void getM(float Y[][3], float IX[][4], /*@out@*/ float M[][3], mwSignedIndex i, mwSignedIndex j, mwSignedIndex k)
/* Determines the affine transform (M) mapping from
   [X+repmat([i j k]', 1,4) ; 1 1 1 1] to [Y ; 1 1 1 1], where IX = inv([X ; 1 1 1 1]);
   This is given by:
        M = Y*inv([X+repmat([i j k]', 1,4) ; 1 1 1 1]);
   or more efficiently by:
        M = Y*(IX - IX(:,1:3)*[i j k]'); */
{
    float ix30, ix31, ix32, ix33;

    ix30 = IX[3][0] - ((float)i*IX[0][0] + (float)j*IX[1][0] + (float)k*IX[2][0]);
    ix31 = IX[3][1] - ((float)i*IX[0][1] + (float)j*IX[1][1] + (float)k*IX[2][1]);
    ix32 = IX[3][2] - ((float)i*IX[0][2] + (float)j*IX[1][2] + (float)k*IX[2][2]);
    ix33 = IX[3][3] - ((float)i*IX[0][3] + (float)j*IX[1][3] + (float)k*IX[2][3]);

    M[0][0] = Y[0][0]*IX[0][0] + Y[1][0]*IX[0][1] + Y[2][0]*IX[0][2] + Y[3][0]*IX[0][3];
    M[0][1] = Y[0][1]*IX[0][0] + Y[1][1]*IX[0][1] + Y[2][1]*IX[0][2] + Y[3][1]*IX[0][3];
    M[0][2] = Y[0][2]*IX[0][0] + Y[1][2]*IX[0][1] + Y[2][2]*IX[0][2] + Y[3][2]*IX[0][3];

    M[1][0] = Y[0][0]*IX[1][0] + Y[1][0]*IX[1][1] + Y[2][0]*IX[1][2] + Y[3][0]*IX[1][3];
    M[1][1] = Y[0][1]*IX[1][0] + Y[1][1]*IX[1][1] + Y[2][1]*IX[1][2] + Y[3][1]*IX[1][3];
    M[1][2] = Y[0][2]*IX[1][0] + Y[1][2]*IX[1][1] + Y[2][2]*IX[1][2] + Y[3][2]*IX[1][3];

    M[2][0] = Y[0][0]*IX[2][0] + Y[1][0]*IX[2][1] + Y[2][0]*IX[2][2] + Y[3][0]*IX[2][3];
    M[2][1] = Y[0][1]*IX[2][0] + Y[1][1]*IX[2][1] + Y[2][1]*IX[2][2] + Y[3][1]*IX[2][3];
    M[2][2] = Y[0][2]*IX[2][0] + Y[1][2]*IX[2][1] + Y[2][2]*IX[2][2] + Y[3][2]*IX[2][3];

    M[3][0] = Y[0][0]*ix30     + Y[1][0]*ix31     + Y[2][0]*ix32     + Y[3][0]*ix33;
    M[3][1] = Y[0][1]*ix30     + Y[1][1]*ix31     + Y[2][1]*ix32     + Y[3][1]*ix33;
    M[3][2] = Y[0][2]*ix30     + Y[1][2]*ix31     + Y[2][2]*ix32     + Y[3][2]*ix33;
}

static void mulMM(float A[][3], float B[][3], float C[][3])
/* [A ; 0 0 0 1] = [B ; 0 0 0 1]*[C ; 0 0 0 1]; */
{
    A[0][0] = B[0][0]*C[0][0] + B[1][0]*C[0][1] + B[2][0]*C[0][2];
    A[0][1] = B[0][1]*C[0][0] + B[1][1]*C[0][1] + B[2][1]*C[0][2];
    A[0][2] = B[0][2]*C[0][0] + B[1][2]*C[0][1] + B[2][2]*C[0][2];

    A[1][0] = B[0][0]*C[1][0] + B[1][0]*C[1][1] + B[2][0]*C[1][2];
    A[1][1] = B[0][1]*C[1][0] + B[1][1]*C[1][1] + B[2][1]*C[1][2];
    A[1][2] = B[0][2]*C[1][0] + B[1][2]*C[1][1] + B[2][2]*C[1][2];

    A[2][0] = B[0][0]*C[2][0] + B[1][0]*C[2][1] + B[2][0]*C[2][2];
    A[2][1] = B[0][1]*C[2][0] + B[1][1]*C[2][1] + B[2][1]*C[2][2];
    A[2][2] = B[0][2]*C[2][0] + B[1][2]*C[2][1] + B[2][2]*C[2][2];

    A[3][0] = B[0][0]*C[3][0] + B[1][0]*C[3][1] + B[2][0]*C[3][2] + B[3][0];
    A[3][1] = B[0][1]*C[3][0] + B[1][1]*C[3][1] + B[2][1]*C[3][2] + B[3][1];
    A[3][2] = B[0][2]*C[3][0] + B[1][2]*C[3][1] + B[2][2]*C[3][2] + B[3][2];
}

static void mulMX(/*@out@*/ float A[][3], float B[][3], float C[][3])
/* [A ; 1 1 1 1] = [B ; 0 0 0 1]*[C ; 1 1 1 1]; */
{
    A[0][0] = B[0][0]*C[0][0] + B[1][0]*C[0][1] + B[2][0]*C[0][2] + B[3][0];
    A[0][1] = B[0][1]*C[0][0] + B[1][1]*C[0][1] + B[2][1]*C[0][2] + B[3][1];
    A[0][2] = B[0][2]*C[0][0] + B[1][2]*C[0][1] + B[2][2]*C[0][2] + B[3][2];

    A[1][0] = B[0][0]*C[1][0] + B[1][0]*C[1][1] + B[2][0]*C[1][2] + B[3][0];
    A[1][1] = B[0][1]*C[1][0] + B[1][1]*C[1][1] + B[2][1]*C[1][2] + B[3][1];
    A[1][2] = B[0][2]*C[1][0] + B[1][2]*C[1][1] + B[2][2]*C[1][2] + B[3][2];

    A[2][0] = B[0][0]*C[2][0] + B[1][0]*C[2][1] + B[2][0]*C[2][2] + B[3][0];
    A[2][1] = B[0][1]*C[2][0] + B[1][1]*C[2][1] + B[2][1]*C[2][2] + B[3][1];
    A[2][2] = B[0][2]*C[2][0] + B[1][2]*C[2][1] + B[2][2]*C[2][2] + B[3][2];

    A[3][0] = B[0][0]*C[3][0] + B[1][0]*C[3][1] + B[2][0]*C[3][2] + B[3][0];
    A[3][1] = B[0][1]*C[3][0] + B[1][1]*C[3][1] + B[2][1]*C[3][2] + B[3][1];
    A[3][2] = B[0][2]*C[3][0] + B[1][2]*C[3][1] + B[2][2]*C[3][2] + B[3][2];
}

static void invertM(float M[][3], /*@out@*/ float IM[][3])
/* [IM ; 0 0 0 1] = inv([M ; 0 0 0 1]); */
{
    float id;
    id = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])+
         M[0][1]*(M[1][2]*M[2][0]-M[1][0]*M[2][2])+
         M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);

    id = 1.0f/id;
    IM[0][0] = (M[1][1]*M[2][2]-M[1][2]*M[2][1])*id;
    IM[0][1] = (M[0][2]*M[2][1]-M[0][1]*M[2][2])*id;
    IM[0][2] = (M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;

    IM[1][0] = (M[1][2]*M[2][0]-M[1][0]*M[2][2])*id;
    IM[1][1] = (M[0][0]*M[2][2]-M[0][2]*M[2][0])*id;
    IM[1][2] = (M[0][2]*M[1][0]-M[0][0]*M[1][2])*id;

    IM[2][0] = (M[1][0]*M[2][1]-M[1][1]*M[2][0])*id;
    IM[2][1] = (M[0][1]*M[2][0]-M[0][0]*M[2][1])*id;
    IM[2][2] = (M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;

    IM[3][0] = (M[1][0]*(M[3][1]*M[2][2]-M[2][1]*M[3][2])+
            M[1][1]*(M[2][0]*M[3][2]-M[3][0]*M[2][2])+
            M[1][2]*(M[3][0]*M[2][1]-M[2][0]*M[3][1]))*id;
    IM[3][1] = (M[0][0]*(M[2][1]*M[3][2]-M[3][1]*M[2][2])+
            M[0][1]*(M[3][0]*M[2][2]-M[2][0]*M[3][2])+
            M[0][2]*(M[2][0]*M[3][1]-M[3][0]*M[2][1]))*id;
    IM[3][2] = (M[0][0]*(M[3][1]*M[1][2]-M[1][1]*M[3][2])+
            M[0][1]*(M[1][0]*M[3][2]-M[3][0]*M[1][2])+
            M[0][2]*(M[3][0]*M[1][1]-M[1][0]*M[3][1]))*id;
}

/****************************************************************************************************/
/* These routines are for locating integer co-ordinates that lie inside a tetrahedron.  See:
    J. Ashburner, J. Andersson and K. J. Friston (2000).
    "Image Registration using a Symmetric Prior - in Three-Dimensions".
    Human Brain Mapping 9(4):212-225. */

static void scan_line(float lin[], mwSignedIndex y, mwSignedIndex z, mwSize *n, /*@out@*/ mwSignedIndex vox[][3], mwSize maxv)
{
    float p[2], t;
    mwSignedIndex x, xe;

    /* sort p into ascending order of x */
    p[0] = lin[0]; p[1] = lin[1];
    if (p[1]<p[0]) {t = p[1]; p[1] = p[0]; p[0] = t;}

    /* find voxels where x is integer */
    for(x=(mwSignedIndex)ceil((double)p[0]), xe=(mwSignedIndex)floor((double)p[1]); x<=xe; x++)
    {
        if ((*n)>=maxv-1)
            mexErrMsgTxt("Too many voxels inside a tetrahedron");
        vox[*n][0] = x;
        vox[*n][1] = y;
        vox[*n][2] = z;
        (*n)++;
    }
}

static void scan_triangle(float tri[][2], mwSignedIndex z, mwSize *n, /*@out@*/ mwSignedIndex vox[][3], mwSize maxv)
{
    float *p[3], *t, lin[2];
    float x1, x2, y1, y2;
    mwSignedIndex y, ye, i;

    /* sort p into ascending order of y */
    p[0] = tri[0]; p[1] = tri[1]; p[2] = tri[2];
    if (p[1][1]<p[0][1]) {t = p[1]; p[1] = p[0]; p[0] = t;}
    if (p[2][1]<p[1][1]) {t = p[2]; p[2] = p[1]; p[1] = t;}
    if (p[1][1]<p[0][1]) {t = p[1]; p[1] = p[0]; p[0] = t;}

    /* find lower lines cutting triangle where y is integer */
    for(y=(mwSignedIndex)ceil((double)p[0][1]), ye=(mwSignedIndex)floor((double)p[1][1]); y<=ye; y++)
    {
        x1 = p[0][0]; y1 = p[0][1];
        for(i=0; i<2; i++)
        {
            x2 = p[i+1][0]; y2 = p[i+1][1];
            if (y2-y1<=0)
                lin[i] = (x1+x2)/2.0f;
            else
                lin[i] = (x1*(y2-(float)y)+x2*((float)y-y1))/(y2-y1);
        }
        scan_line(lin,y,z, n,vox,maxv);
    }

    /* find upper lines cutting triangle where y is integer */
    for(y=(mwSignedIndex)ceil((double)p[1][1]), ye=(mwSignedIndex)floor((double)p[2][1]); y<=ye; y++)
    {
        x2 = p[2][0]; y2 = p[2][1];
        for(i=0; i<2; i++)
        {
            x1 = p[i][0]; y1 = p[i][1];
            if (y2-y1<=0)
                lin[i] = (x1+x2)/2.0f;
            else
                lin[i] = (x1*(y2-(float)y)+x2*((float)y-y1))/(y2-y1);
        }
        scan_line(lin,y,z, n,vox,maxv);
    }
}


static void scan_tetrahedron(float Y[][3], /*@out@*/ mwSize *n, /*@out@*/ mwSignedIndex vox[][3], mwSize maxv)
/* Y are the vertices of the tetrahedron.  n are the number of located
   integer co-ordinates, vox are the co-ordinates found and maxv are the
   maximum number of co-ordinates allowed. */
{
    float *p[4], *t, tri[4][2], x1, x2, y1, y2, z1, z2;
    mwSignedIndex z, ze, i;

    *n = 0;

    /* sort p into ascending order of z */
    p[0] = Y[0]; p[1] = Y[1]; p[2] = Y[2]; p[3] = Y[3];
    if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}
    if (p[2][2]<p[1][2]) {t = p[2]; p[2] = p[1]; p[1] = t;}
    if (p[3][2]<p[2][2]) {t = p[3]; p[3] = p[2]; p[2] = t;}
    if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}
    if (p[2][2]<p[1][2]) {t = p[2]; p[2] = p[1]; p[1] = t;}
    if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}

    /* find lower triangles that intersect tetrahedron where z is integer */
    for(z=(mwSignedIndex)ceil((double)p[0][2]), ze=(mwSignedIndex)floor((double)p[1][2]); z<=ze; z++)
    {
        x1 = p[0][0]; y1 = p[0][1]; z1 = p[0][2];
        for(i=0; i<3; i++)
        {
            x2 = p[i+1][0]; y2 = p[i+1][1]; z2 = p[i+1][2];
            if (z2-z1<=0)
            {
                tri[i][0] = (x1+x2)/2.0f;
                tri[i][1] = (y1+y2)/2.0f;
            }
            else
            {
                float  t2 = z2-(float)z, t1 = (float)z-z1, tmp = z2-z1;
                tri[i][0] = (x1*t2+x2*t1)/tmp;
                tri[i][1] = (y1*t2+y2*t1)/tmp;
            }
        }
        scan_triangle(tri,z, n,vox,maxv);
    }

    /* find quadrilaterals that intersect tetrahedron where z is integer */
    /* each quadrilateral divided into two triangles */
    for(z=(mwSignedIndex)ceil((double)p[1][2]), ze=(mwSignedIndex)floor((double)p[2][2]); z<=ze; z++)
    {
        static int ii[] = {0,1,1,0}, jj[] = {3,3,2,2};

        for(i=0; i<4; i++)
        {
            x1 = p[ii[i]][0]; y1 = p[ii[i]][1]; z1 = p[ii[i]][2];
            x2 = p[jj[i]][0]; y2 = p[jj[i]][1]; z2 = p[jj[i]][2];
            if (z2-z1<=0)
            {
                tri[i][0] = (x1+x2)/2.0f;
                tri[i][1] = (y1+y2)/2.0f;
            }
            else
            {
                float  t2 = z2-(float)z, t1 = (float)z-z1, tmp = z2-z1;
                tri[i][0] = (x1*t2+x2*t1)/tmp;
                tri[i][1] = (y1*t2+y2*t1)/tmp;
            }
        }
        scan_triangle(tri,z, n,vox,maxv);
        tri[1][0] = tri[3][0];
        tri[1][1] = tri[3][1];
        scan_triangle(tri,z, n,vox,maxv);
    }

    /* find upper triangles that intersect tetrahedron where z is integer */
    for(z=(mwSignedIndex)ceil((double)p[2][2]), ze=(mwSignedIndex)floor((double)p[3][2]); z<=ze; z++)
    {
        x2 = p[3][0]; y2 = p[3][1]; z2 = p[3][2];
        for(i=0; i<3; i++)
        {
            x1 = p[i][0]; y1 = p[i][1]; z1 = p[i][2];
            if (z2-z1<=0)
            {
                tri[i][0] = (x1+x2)/2.0f;
                tri[i][1] = (y1+y2)/2.0f;
            }
            else
            {
                float t2 = z2-(float)z, t1 = (float)z-z1, tmp = z2-z1;
                tri[i][0] = (x1*t2+x2*t1)/tmp;
                tri[i][1] = (y1*t2+y2*t1)/tmp;
            }
        }
        scan_triangle(tri,z, n,vox,maxv);
    }
}

/****************************************************************************************************/

/* Division of a cube into two alternating patterns of tetrahedra.
  This pattern is repeated in a 3D checkerboard pattern, such that
  the whole volume is covered. */
static float x[2][5][4][3] = {
{{{ 0,0,0},{ 1,0,1},{ 1,0,0},{ 1,1,0}},
 {{ 0,0,0},{ 1,0,1},{ 0,1,1},{ 0,0,1}},
 {{ 0,0,0},{ 0,1,0},{ 0,1,1},{ 1,1,0}},
 {{ 0,0,0},{ 1,0,1},{ 1,1,0},{ 0,1,1}},
 {{ 1,1,1},{ 1,1,0},{ 0,1,1},{ 1,0,1}},},

{{{ 1,0,0},{ 0,0,1},{ 0,0,0},{ 0,1,0}},
 {{ 1,0,0},{ 0,0,1},{ 1,1,1},{ 1,0,1}},
 {{ 1,0,0},{ 1,1,0},{ 1,1,1},{ 0,1,0}},
 {{ 1,0,0},{ 0,0,1},{ 0,1,0},{ 1,1,1}},
 {{ 0,1,1},{ 0,1,0},{ 1,1,1},{ 0,0,1}},},
};


/* Set up matrices (IX) for each tetrahedron, so that future computations can be
   made faster.  Also set up a few relative file offsets. */
static float ix[2][5][4][4];
static mwSignedIndex off[2][4][5];
static void setup_consts(mwSize dim[])
{
    mwSignedIndex i, j, k;
    for(k=0; k<2; k++)
        for(i=0; i<5; i++)
        {
            invertX(x[k][i], ix[k][i]);
            for(j=0; j<4; j++)
                off[k][j][i] = (mwSignedIndex)(x[k][i][j][0]+dim[0]*(x[k][i][j][1]+dim[1]*x[k][i][j][2]));
        }
}

/* Compute the inverse deformation field within a single cube */
static void invert_it(mwSignedIndex x0, mwSignedIndex x1, mwSignedIndex x2, float *y0, float *y1, float *y2,
    mwSize dim_iy[], float *iy0, float *iy1, float *iy2, float M1[][3], float M2[][3], int pass)
{
    mwSignedIndex i, j, vox[MAXV][3];
    float  Y0[4][3], Y[4][3], M[4][3], IM[4][3];
    mwSize nvox;
    int k;

    /* Determine tetrahedral arrangement */
    k = (int)((x0%2)==(mwSignedIndex)((x1%2)==(x2%2)));
    if (pass==1) k = (int)(k==0);
 
    for(i=0; i<5; i++) /* Five tetrahedra within a cube */
    {
        /* Find the vertices (in mm space) */
        Y0[0][0] = y0[off[k][0][i]]; Y0[0][1] = y1[off[k][0][i]]; Y0[0][2] = y2[off[k][0][i]];
        Y0[1][0] = y0[off[k][1][i]]; Y0[1][1] = y1[off[k][1][i]]; Y0[1][2] = y2[off[k][1][i]];
        Y0[2][0] = y0[off[k][2][i]]; Y0[2][1] = y1[off[k][2][i]]; Y0[2][2] = y2[off[k][2][i]];
        Y0[3][0] = y0[off[k][3][i]]; Y0[3][1] = y1[off[k][3][i]]; Y0[3][2] = y2[off[k][3][i]];

        /* Convert vertex co-ordinates to voxels */
        mulMX(Y, M1, Y0);

        /* Compute affine transform mapping vertices together */
        getM(Y, ix[k][i], M, x0, x1, x2);

        if (mxIsFinite(M[0][0])) /* Prevent from bombing out when NaNs are encountered */
        {
            /* Find integer co-ordinates within tetrahedron */
            scan_tetrahedron(Y, &nvox, vox, MAXV);

            if (nvox>0)
            {
                /* Invert the affine mapping */
                invertM(M, IM);

                /* Convert the mapping from voxels to mm */
                mulMM(M, M2, IM);

                if (pass==0)
                {
                    /* Insert the new mappings into each voxel within the tetrahedron */
                    for(j=0; j<(mwSignedIndex)nvox; j++)
                    {
                        if ((vox[j][0]>=1) && (vox[j][0]<=(mwSignedIndex)dim_iy[0]) &&
                            (vox[j][1]>=1) && (vox[j][1]<=(mwSignedIndex)dim_iy[1]) &&
                            (vox[j][2]>=1) && (vox[j][2]<=(mwSignedIndex)dim_iy[2]))
                        {
                            mwSignedIndex o  = vox[j][0]+dim_iy[0]*(vox[j][1]+dim_iy[1]*vox[j][2]);

                            iy0[o] = M[0][0]*(float)vox[j][0] + M[1][0]*(float)vox[j][1] + M[2][0]*(float)vox[j][2] + M[3][0];
                            iy1[o] = M[0][1]*(float)vox[j][0] + M[1][1]*(float)vox[j][1] + M[2][1]*(float)vox[j][2] + M[3][1];
                            iy2[o] = M[0][2]*(float)vox[j][0] + M[1][2]*(float)vox[j][1] + M[2][2]*(float)vox[j][2] + M[3][2];
                        }
                    }
                }
                else
                {
                    /* Average the new mappings with those from the 0th pass */
                    for(j=0; j<(mwSignedIndex)nvox; j++)
                    {
                        if ((vox[j][0]>=1) && (vox[j][0]<=(mwSignedIndex)dim_iy[0]) &&
                            (vox[j][1]>=1) && (vox[j][1]<=(mwSignedIndex)dim_iy[1]) &&
                            (vox[j][2]>=1) && (vox[j][2]<=(mwSignedIndex)dim_iy[2]))
                        {
                            mwSignedIndex o  = vox[j][0]+dim_iy[0]*(vox[j][1]+dim_iy[1]*vox[j][2]);

                            iy0[o] = (iy0[o] + M[0][0]*(float)vox[j][0] + M[1][0]*(float)vox[j][1] + M[2][0]*(float)vox[j][2] + M[3][0])/2.0f;
                            iy1[o] = (iy1[o] + M[0][1]*(float)vox[j][0] + M[1][1]*(float)vox[j][1] + M[2][1]*(float)vox[j][2] + M[3][1])/2.0f;
                            iy2[o] = (iy2[o] + M[0][2]*(float)vox[j][0] + M[1][2]*(float)vox[j][1] + M[2][2]*(float)vox[j][2] + M[3][2])/2.0f;
                        }
                    }
                }
            }
        }
    }
}


/* Some regions of the inverse deformation may be undefined. */
static void setnan(float *dat, mwSize n)
{
    mwIndex j;
    float NaN;
    NaN = mxGetNaN();
    for (j=0; j<n; j++) dat[j] = NaN;
}

void invdef(mwSize *dim_y,  float  y0[],
                   mwSize *dim_iy, float iy0[],  float M1[][3],  float M2[][3])
{
    mwSignedIndex x2, x1, x0;
    int pass;
    float *y1=0, *y2=0, *iy1=0, *iy2=0;

    setup_consts(dim_y);
    setnan(iy0, dim_iy[0]*dim_iy[1]*dim_iy[2]*3);

    /* Convert arrays such that the first voxel is at [1 1 1] rather than [0 0 0]
       (a trick used by f2c). Generate pointers to second and third componsnts. */
    y0  -= 1+dim_y[0]*(1 + dim_y[1]);
    y1   = y0+dim_y[0]*dim_y[1]*dim_y[2];
    y2   = y0+dim_y[0]*dim_y[1]*dim_y[2]*2;

    iy0 -= 1+dim_iy[0]*(1 + dim_iy[1]);
    iy1  = iy0+dim_iy[0]*dim_iy[1]*dim_iy[2];
    iy2  = iy0+dim_iy[0]*dim_iy[1]*dim_iy[2]*2;

    /* Two passes because there are two possible tetrahedral arrangements */
    for (pass=0; pass<=1; pass++)
    {
        /* Loop over all cubes in the deformation field. */
        for(x2=1; x2<(mwSignedIndex)dim_y[2]; x2++)
        {
            for(x1=1; x1<(mwSignedIndex)dim_y[1]; x1++)
            {
                for(x0=1; x0<(mwSignedIndex)dim_y[0]; x0++)
                {
                    mwSignedIndex o = x0 + dim_y[0]*(x1 + x2*dim_y[1]);
                    invert_it(x0, x1, x2, y0+o, y1+o, y2+o, dim_iy, iy0, iy1, iy2, M1, M2, pass);
                }
            }
        }
    }
}

