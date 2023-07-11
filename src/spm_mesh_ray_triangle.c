/*
 * Guillaume Flandin
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include "mex.h"
#include "external/raytri/raytri.c"

/* Gateway Function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *tri = NULL, *ray = NULL;
    double *H = NULL, *P = NULL;
    mwIndex *h = NULL;
    double t, u, v;
    mwSize nv, nr;
    mwIndex i, j, n = 32;
    int hit;
    double v1[3], v2[3], v3[3];
    double orig[3], dir[3];

    if (nrhs < 2) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 2) mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 2) mexErrMsgTxt("Too many output arguments.");

    if (!mxIsDouble(prhs[0]) || mxGetM(prhs[0]) != 9)
        mexErrMsgTxt("Triangles have to be provided as a 9xNV array.");
    nv  = mxGetN(prhs[0]);
    tri = mxGetPr(prhs[0]);

    if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 6)
        mexErrMsgTxt("Rays have to be provided as a 6xNR array.");
    nr = mxGetN(prhs[1]);
    ray = mxGetPr(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(nr, n, mxREAL);
    H = mxGetPr(plhs[0]);
    h = mxCalloc(nr, sizeof(mwIndex));
    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(nr, n, mxREAL);
        P = mxGetPr(plhs[1]);
        /* should initialise to -Inf */
    }

    for (i=0;i<nv;i++)
    {
        v1[0] = tri[0+9*i]; v1[1] = tri[1+9*i]; v1[2] = tri[2+9*i];
        v2[0] = tri[3+9*i]; v2[1] = tri[4+9*i]; v2[2] = tri[5+9*i];
        v3[0] = tri[6+9*i]; v3[1] = tri[7+9*i]; v3[2] = tri[8+9*i];

        for (j=0;j<nr;j++)
        {
            orig[0] = ray[0+6*j]; orig[1] = ray[1+6*j]; orig[2] = ray[2+6*j];
            dir[0]  = ray[3+6*j]; dir[1]  = ray[4+6*j]; dir[2]  = ray[5+6*j];

            hit = intersect_triangle1(orig, dir, v1, v2, v3, &t, &u, &v);

            if (hit) /* could store positive hits only */
            {
                H[j+nr*h[j]] = i + 1;
                if (nlhs > 1) P[j+nr*h[j]] = t;
                h[j] += 1;
                if (h[j] == n)
                {
                    mexWarnMsgTxt("Max number of intersections reached. Reallocating.");
                    n *= 2;
                    H = mxRealloc(H, nr*n*sizeof(double));
                    mxSetPr(plhs[0],H); mxSetN(plhs[0],n);
                    if (nlhs > 1) 
                    {
                        P = mxRealloc(P, nr*n*sizeof(double));
                        mxSetPr(plhs[1],P); mxSetN(plhs[1],n);
                    }
                }
            }
        }
    }
    mxFree(h);
}
