/*
 * $Id: spm_resels_vol.c 4453 2011-09-02 10:47:25Z guillaume $
 * John Ashburner
 */
 
/*
See Worsley et al. (1996), Human Brain Mapping 4:58-73 for a description
of what it does.
*/

#include <math.h>
#include "mex.h"
#include "spm_mapping.h"

static void resel_fun(int *curr, int *prev, /* current and previous planes */
    int m, int n, /* image dimensions */
    int *P,   /* # points */
    int E[3], /* # edges  */
    int F[4], /* # faces  */
    int *C)   /* # cubes  */
{
    int p=0,ex=0,ey=0,ez=0,fxy=0,fxz=0,fyz=0,c=0;
    int i, j;
    for(i=1; i<m; i++)
        for(j=1; j<n; j++)
        {
            int o = i+j*m;

            if (curr[o])
            {
                p ++;

                /* A simple way of doing it... 
                if (curr[o-1]) ex ++;
                if (curr[o-m]) ey ++;
                if (prev[o  ]) ez ++;
                if (curr[o-1] && curr[o-m] && curr[o-m-1]) fxy ++;
                if (curr[o-1] && prev[o  ] && prev[o  -1]) fxz ++;
                if (curr[o-m] && prev[o  ] && prev[o  -m]) fyz ++;
                if (curr[o-1] && curr[o-m] && curr[o-m-1] && prev[o]
                 && prev[o-1] && prev[o-m] && prev[o-m-1]) c ++;
                */

                /* It could be optimized much further with "if else" statements
                   but it was beginning to hurt my head */
                if (curr[o-1])
                {
                    ex ++;
                    if (curr[o-m] && curr[o-m-1])
                    {
                        fxy ++;
                        if (prev[o] && prev[o-1] && prev[o-m] && prev[o-m-1]) c ++;
                    }
                }
                if (curr[o-m])
                {
                    ey ++;
                    if (prev[o  ] && prev[o  -m]) fyz ++;
                }
                if (prev[o  ])
                {
                    ez ++;
                    if (curr[o-1] && prev[o  -1]) fxz ++;
                }
            }
        }
    *P   += p;
    E[0] += ex;
    E[1] += ey;
    E[2] += ez;
    F[0] += fxy;
    F[1] += fxz;
    F[2] += fyz;
    *C   += c;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    MAPTYPE *map, *get_maps();
    mwSize m,n,k;
    mwIndex i;
    int nn, E[3], F[3], P, C;
    double *R, r[3], *img;
    int *curr, *prev, *tmpp;

    if (nrhs != 2 || nlhs > 1)
    {
        mexErrMsgTxt("Incorrect usage.");
    }

    map = get_maps(prhs[0], &nn);
    if (nn!=1)
    {
        free_maps(map, nn);
        mexErrMsgTxt("Bad image handle dimensions.");
    }

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || !mxIsDouble(prhs[1]))
    {
        free_maps(map, 1);
        mexErrMsgTxt("Second argument must be numeric, real, full and double.");
    }

    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 3)
    {
        free_maps(map, nn);
        mexErrMsgTxt("Second argument must contain three elements.");
    }

    r[0] = 1.0/mxGetPr(prhs[1])[0];
    r[1] = 1.0/mxGetPr(prhs[1])[1];
    r[2] = 1.0/mxGetPr(prhs[1])[2];

    plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);
    R       = mxGetPr(plhs[0]);

    m = map->dim[0];
    n = map->dim[1];
    k = map->dim[2];

    curr = (int *)mxCalloc((m+1)*(n+1),sizeof(int)); /* current plane */
    prev = (int *)mxCalloc((m+1)*(n+1),sizeof(int)); /* previous plane */
    img  = (double *)mxCalloc((m+1)*(n+1),sizeof(double));

    P = C = E[0] = E[1] = E[2] = F[0] = F[1] = F[2] = 0;
    for(i=0; i<k; i++)
    {
        int i1, j1;
        static double mat[16] = {
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0};
        mat[14] = i+1.0;
        slice(mat, img, m, n, map, 0, 0);
        for(i1=0;i1<m; i1++)
            for(j1=0; j1<n; j1++)
                if (mxIsFinite(img[i1+m*j1]) && img[i1+m*j1])
                    curr[i1+1+(m+1)*(j1+1)] = 1;
                else
                    curr[i1+1+(m+1)*(j1+1)] = 0;

        /* count edges, faces etc */
        resel_fun(curr,prev,m+1,n+1, &P, E, F, &C);

        /* make current plane previous */
        tmpp = prev; prev = curr; curr = tmpp;
    }

    (void)mxFree((char *)curr);
    (void)mxFree((char *)prev);
    (void)mxFree((char *)img);
    free_maps(map, 1);

    R[0] = P - (E[0]+E[1]+E[2])+(F[0]+F[1]+F[2])-C;
    R[1] = (E[0]-F[0]-F[1]+C)*r[0] + (E[1]-F[0]-F[2]+C)*r[1] + (E[2]-F[1]-F[2]+C)*r[2];
    R[2] = (F[0]-C)*r[0]*r[1] + (F[1]-C)*r[0]*r[2] + (F[2]-C)*r[1]*r[2];
    R[3] = C*r[0]*r[1]*r[2];
}
