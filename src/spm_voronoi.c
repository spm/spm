/*
 * $Id: spm_voronoi.c 2052 2008-09-08 16:30:57Z guillaume $
 * Guillaume Flandin
 */

#include <math.h>
#include <string.h>
#include "mex.h"

/* 3D Point Structure */
typedef struct {
    int x;
    int y;
    int z;
} point3d;

/* Chamfer Mask Structure */
typedef struct {
    int x;
    int y;
    int z;
    double v;
} mask;

/* Bucket Structure */
typedef struct {
    point3d *pt;
    int size;
    int memsize;
} bucket;

/* Queue Structure */
typedef struct {
    bucket *b;
    int nbbucket;
} queue;

/* d4 mask definition */
mask maskd4[] = {{ 0,  0, -1,  1.0},
                 { 0, -1,  0,  1.0},
                 {-1,  0,  0,  1.0},
                 { 1,  0,  0,  1.0},
                 { 0,  1,  0,  1.0},
                 { 0,  0,  1,  1.0}};

/* d8 mask definition */
mask maskd8[] = {{-1, -1,  0,  1.0},
                 {-1,  0, -1,  1.0},
                 {-1,  0,  0,  1.0},
                 {-1,  0,  1,  1.0},
                 {-1,  1,  0,  1.0},
                 { 0, -1, -1,  1.0},
                 { 0, -1,  0,  1.0},
                 { 0,  0, -1,  1.0},
                 { 0, -1,  1,  1.0},
                 { 1,  1,  0,  1.0},
                 { 1,  0,  1,  1.0},
                 { 1,  0,  0,  1.0},
                 { 1,  0, -1,  1.0},
                 { 1, -1,  0,  1.0},
                 { 0,  1,  1,  1.0},
                 { 0,  1,  0,  1.0},
                 { 0,  0,  1,  1.0},
                 { 0,  1, -1,  1.0},
                 };
    
/* d333 definition : 3 4 5*/
mask maskd34[] = {{-1, -1, -1,  5.0},
                  {-1, -1,  0,  4.0},
                  {-1, -1,  1,  5.0},
                  {-1,  0, -1,  4.0},
                  {-1,  0,  0,  3.0},
                  {-1,  0,  1,  4.0},
                  {-1,  1, -1,  5.0},
                  {-1,  1,  0,  4.0},
                  {-1,  1,  1,  5.0},
                  { 0, -1, -1,  4.0},
                  { 0, -1,  0,  3.0},
                  { 0,  0, -1,  3.0},
                  { 0, -1,  1,  4.0},
                  { 1,  1,  1,  5.0},
                  { 1,  1,  0,  4.0},
                  { 1,  1, -1,  5.0},
                  { 1,  0,  1,  4.0},
                  { 1,  0,  0,  3.0},
                  { 1,  0, -1,  4.0},
                  { 1, -1,  1,  5.0},
                  { 1, -1,  0,  4.0},
                  { 1, -1, -1,  5.0},
                  { 0,  1,  1,  4.0},
                  { 0,  1,  0,  3.0},
                  { 0,  0,  1,  3.0},
                  { 0,  1, -1,  4.0},
                  };

/* d5711 mask definition : 5 7 9 11 12 15  */
mask maskd5711[] = {{-1,  -2,  -2,  15.0},
                    { 1,  -2,  -2,  15.0},
                    {-2,  -1,  -2,  15.0},
                    {-1,  -1,  -2,  12.0},
                    { 0,  -1,  -2,  11.0},
                    { 1,  -1,  -2,  12.0},
                    { 2,  -1,  -2,  15.0},
                    {-1,   0,  -2,  11.0},
                    { 1,   0,  -2,  11.0},
                    {-2,   1,  -2,  15.0},
                    {-1,   1,  -2,  12.0},
                    { 0,   1,  -2,  11.0},
                    { 1,   1,  -2,  12.0},
                    { 2,   1,  -2,  15.0},
                    {-1,   2,  -2,  15.0},
                    { 1,   2,  -2,  15.0},
                    
                    {-2,  -2,  -1,  15.0},
                    {-1,  -2,  -1,  12.0},
                    { 0,  -2,  -1,  11.0},
                    { 1,  -2,  -1,  12.0},
                    { 2,  -2,  -1,  15.0},
                    {-2,  -1,  -1,  12.0},
                    {-1,  -1,  -1,   9.0},
                    { 0,  -1,  -1,   7.0},
                    { 1,  -1,  -1,   9.0},
                    { 2,  -1,  -1,  12.0},
                    {-2,   0,  -1,  11.0},
                    {-1,   0,  -1,   7.0},
                    { 0,   0,  -1,   5.0},
                    { 1,   0,  -1,   7.0},
                    { 2,   0,  -1,  11.0},
                    {-2,   1,  -1,  12.0},
                    {-1,   1,  -1,   9.0},
                    { 0,   1,  -1,   7.0},
                    { 1,   1,  -1,   9.0},
                    { 2,   1,  -1,  12.0},
                    {-2,   2,  -1,  15.0},
                    {-1,   2,  -1,  12.0},
                    { 0,   2,  -1,  11.0},
                    { 1,   2,  -1,  12.0},
                    { 2,   2,  -1,  15.0},
                    
                    {-1,  -2,   0,  11.0},
                    { 1,  -2,   0,  11.0},
                    {-2,  -1,   0,  11.0},
                    {-1,  -1,   0,   7.0},
                    { 0,  -1,   0,   5.0},
                    { 1,  -1,   0,   7.0},
                    { 2,  -1,   0,  11.0},
                    {-1,   0,   0,   5.0},
                    { 1,   0,   0,   5.0},
                    {-2,   1,   0,  11.0},
                    {-1,   1,   0,   7.0},
                    { 0,   1,   0,   5.0},
                    { 1,   1,   0,   7.0},
                    { 2,   1,   0,  11.0},
                    {-1,   2,   0,  11.0},
                    { 1,   2,   0,  11.0},
                    
                    {-2,  -2,   1,  15.0},
                    {-1,  -2,   1,  12.0},
                    { 0,  -2,   1,  11.0},
                    { 1,  -2,   1,  12.0},
                    { 2,  -2,   1,  15.0},
                    {-2,  -1,   1,  12.0},
                    {-1,  -1,   1,   9.0},
                    { 0,  -1,   1,   7.0},
                    { 1,  -1,   1,   9.0},
                    { 2,  -1,   1,  12.0},
                    {-2,   0,   1,  11.0},
                    {-1,   0,   1,   7.0},
                    { 0,   0,   1,   5.0},
                    { 1,   0,   1,   7.0},
                    { 2,   0,   1,  11.0},
                    {-2,   1,   1,  12.0},
                    {-1,   1,   1,   9.0},
                    { 0,   1,   1,   7.0},
                    { 1,   1,   1,   9.0},
                    { 2,   1,   1,  12.0},
                    {-2,   2,   1,  15.0},
                    {-1,   2,   1,  12.0},
                    { 0,   2,   1,  11.0},
                    { 1,   2,   1,  12.0},
                    { 2,   2,   1,  15.0},
                    
                    {-1,  -2,   2,  15.0},
                    { 1,  -2,   2,  15.0},
                    {-2,  -1,   2,  15.0},
                    {-1,  -1,   2,  12.0},
                    { 0,  -1,   2,  11.0},
                    { 1,  -1,   2,  12.0},
                    { 2,  -1,   2,  15.0},
                    {-1,   0,   2,  11.0},
                    { 1,   0,   2,  11.0},
                    {-2,   1,   2,  15.0},
                    {-1,   1,   2,  12.0},
                    { 0,   1,   2,  11.0},
                    { 1,   1,   2,  12.0},
                    { 2,   1,   2,  15.0},
                    {-1,   2,   2,  15.0},
                    { 1,   2,   2,  15.0}};

/* --- GEODESIC VORONOI SUBROUTINE --- */
void geodesic_voronoi(double *img,         /* domain image mask */
                      unsigned int *size,  /* size of domain image */
                      double *voxsize,     /* size of voxels (in mm) */
                      double *seeds,       /* positions of the seeds */
                      int nbseeds,         /* number of seeds*/
                      double *vor,         /* voronoi diagram */
                      double *dmap,        /* geodesic distance map */
                      char *dist) {        /* type of distance */
    int i, j, a;
    int k = 1, l = 1, p = 1, s = 1;
    unsigned int m = size[0], n = size[1], o = size[2];
    
    mask *cmask = NULL;
    queue q;
    int notempty = 1;
    int current_buck = 0;
    int x, y, z;
    double newd;
    
    /* --- CHAMFER MASK INITIALIZATION --- */
    if (!strcmp(dist,"d4")) {
        cmask = maskd4;
        k = l = p = 1;
        s = 6;
    }
    else if (!strcmp(dist,"d8")) {
        cmask = maskd8;
        k = l = p = 1;
        s = 18;
    }
    else if (!strcmp(dist,"d34")) {
        cmask = maskd34;
        k = l = p = 1;
        s = 26;
    }
    else if (!strcmp(dist,"d5711")) {
        cmask = maskd5711;
        k = l = p = 2;
        s = 100;
        /* mexWarnMsgTxt("Geodesic distances may be underestimated with d5711."); */
    }
    else {
        mexErrMsgTxt("Unknown distance type.");
    }
    
    /* --- VOXEL SIZE WEIGHTING --- */
    
    
    /* --- BUCKET, DISTMAP AND VORONOI INITIALIZATION --- */
    q.b = mxMalloc(sizeof(bucket));
    q.nbbucket = 1;
    
    /* distmap and voronoi initialization */
    for (i=0;i<m*n*o;i++)
        if (img[i] > 0.0) {
            dmap[i] = mxGetInf();  /* inside domain */
            vor[i]  = mxGetNaN();
        }
        else {
            dmap[i] = -1; /* outside domain */
            vor[i]  = -1;
        }

    /* bucket initialization and distmap and voronoi update */
    q.b[0].size = nbseeds;
    q.b[0].pt = mxMalloc(nbseeds*sizeof(point3d));
    q.b[0].memsize = nbseeds;
    for (i=0;i<nbseeds;i++) {
        q.b[0].pt[i].x = (int)(seeds[i]+0.5) - 1;
        q.b[0].pt[i].y = (int)(seeds[i + nbseeds]+0.5) - 1;
        q.b[0].pt[i].z = (int)(seeds[i + 2 * nbseeds]+0.5) - 1;
        if ((q.b[0].pt[i].x < 0) || (q.b[0].pt[i].x >= m) ||
            (q.b[0].pt[i].y < 0) || (q.b[0].pt[i].y >= n) ||
            (q.b[0].pt[i].z < 0) || (q.b[0].pt[i].z >= o))
            mexErrMsgTxt("A seed appears to be outside image.");
        dmap[q.b[0].pt[i].x+m*q.b[0].pt[i].y+m*n*q.b[0].pt[i].z] = 0.0;
        vor[q.b[0].pt[i].x+m*q.b[0].pt[i].y+m*n*q.b[0].pt[i].z]  = i+1;
    }

    /* --- BUCKET SORTING --- */
    while (notempty) {

        /* loop over the not empty bucket with the smallest value */
        for (i=0;i<q.b[current_buck].size;i++) {
            
            /* loop over the chamfer mask */
            for (a=0;a<s;a++) {
                x = q.b[current_buck].pt[i].x + cmask[a].x;
                y = q.b[current_buck].pt[i].y + cmask[a].y;
                z = q.b[current_buck].pt[i].z + cmask[a].z;
                
                /* point inside image ? */
                if ((x >= 0) && (y >= 0) && (z >= 0)
                     && (x < m) && (y < n) && (z < o)) {

                    /* point inside domain ? */
                    if (dmap[x+m*y+m*n*z] > 0.0) {

                        newd = dmap[q.b[current_buck].pt[i].x
                                + m*q.b[current_buck].pt[i].y
                                + m*n*q.b[current_buck].pt[i].z]
                                + cmask[a].v;

                        /* shorter distance ? */
                        if (newd < dmap[x+m*y+m*n*z]) {
                            
                            /* new value for that distance */
                            dmap[x+m*y+m*n*z] = newd;
                            
                            /* Update Voronoi diagram */
                            vor[x+m*y+m*n*z] = vor[q.b[current_buck].pt[i].x
                                    + m*q.b[current_buck].pt[i].y
                                    + m*n*q.b[current_buck].pt[i].z];
                            
                            /* add in bucket corresponding to newd */
                            /* queue size update */
                            if (q.nbbucket <= (int)newd) {
                                q.b = mxRealloc(q.b,((int)newd+1) * sizeof(bucket));
                                for (j=q.nbbucket;j<(int)newd+1;j++) {
                                    q.b[j].size = 0;
                                    q.b[j].memsize = 2;
                                    q.b[j].pt = mxMalloc(2*sizeof(point3d));
                                }
                                q.nbbucket = (int)newd+1;
                            }
                            
                            /* bucket size update */
                            if (q.b[(int)newd].size >= q.b[(int)newd].memsize) {
                                q.b[(int)newd].pt = mxRealloc(q.b[(int)newd].pt,q.b[(int)newd].memsize*2*sizeof(point3d));
                                q.b[(int)newd].memsize *=2;
                            }   
                            q.b[(int)newd].pt[q.b[(int)newd].size].x = x;
                            q.b[(int)newd].pt[q.b[(int)newd].size].y = y;
                            q.b[(int)newd].pt[q.b[(int)newd].size].z = z;
                            q.b[(int)newd].size++;
                        }
                    }
                }
            }
        }
        mxFree(q.b[current_buck].pt);
        /* empty buckets ? */
        if (q.nbbucket == current_buck + 1)
            notempty = 0;
        else
            current_buck++;
    }
    return;
}

/* --- GATEWAY FUNCTION --- */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

double *img = NULL, *distmap = NULL, *vor = NULL;
double *voxsize = NULL, *seeds = NULL;
unsigned int size[3], sdist, nbseeds;
char *dist = NULL;
mxArray *dmap = NULL, *vord = NULL;

/* Check for proper number of arguments. */
if (nrhs < 2)
    mexErrMsgTxt("Not enough input arguments.");
else if (nrhs > 3)
    mexErrMsgTxt("Too many input arguments.");
else if (nlhs > 2)
    mexErrMsgTxt("Too many output arguments.");

/* The input IMG must be a real double array. */
if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    !(mxGetNumberOfDimensions(prhs[0])==2 || 
    mxGetNumberOfDimensions(prhs[0])==3))
    mexErrMsgTxt("Input img must be a real double array.");
size[0] = mxGetDimensions(prhs[0])[0];
size[1] = mxGetDimensions(prhs[0])[1];
if (mxGetNumberOfDimensions(prhs[0])==3)
    size[2] = mxGetDimensions(prhs[0])[2];
else 
    size[2] = 1;
img = mxGetPr(prhs[0]);

/* The input SEEDS must be a real double array [n x 3] */
if (!mxIsDouble(prhs[1]) || 
     mxIsComplex(prhs[1]) || 
     mxGetNumberOfDimensions(prhs[1]) != 2 ||
     mxGetDimensions(prhs[1])[1] != 3)
    mexErrMsgTxt("Input seeds must be a real double array [n x 3].");
seeds = mxGetPr(prhs[1]);
nbseeds = mxGetDimensions(prhs[1])[0];

/* The input DIST must be a string */
if (nrhs > 2) {
    if (!mxIsChar(prhs[2]))
        mexErrMsgTxt("Input dist must be a string.");
    sdist = (mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
    dist = mxCalloc(sdist, sizeof(char));
    mxGetString(prhs[2], dist, sdist);
}
else {
    dist = "d34";
}

/* The input VOXSIZE must be a [x y z] array */
if (nrhs > 3) {
    if (!mxIsDouble(prhs[3]) 
            || mxIsComplex(prhs[3])
            || mxGetNumberOfElements(prhs[3])!=3)
    mexErrMsgTxt("Input voxsize must be a 3 elements vector.");
    voxsize = mxGetPr(prhs[3]);
}
else {
    voxsize = mxMalloc(3*sizeof(double));
    voxsize[0] = voxsize[1] = voxsize[2] = 1.0;
}

/* Create mxArray's */
dmap = mxCreateNumericArray(3, (const int*)size, mxDOUBLE_CLASS, mxREAL);
distmap = mxGetPr(dmap);

vord = mxCreateNumericArray(3, (const int*)size, mxDOUBLE_CLASS, mxREAL);
vor = mxGetPr(vord);

/* Call the Geodesic Voronoi subroutine. */
geodesic_voronoi(img, size, voxsize, seeds, nbseeds, vor, distmap, dist);

/* Assign pointers to output */
if (nlhs > 0) plhs[0] = vord;
else mxDestroyArray(vord);

if (nlhs > 1) plhs[1] = dmap;
else mxDestroyArray(dmap);

/* Free memory */
if (nrhs < 4) mxFree(voxsize);

}
