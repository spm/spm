/*
 * $Id: spm_global.c 4921 2012-09-13 11:16:21Z guillaume $
 * John Ashburner
 */

#include "mex.h"
#include "spm_mapping.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize i, j, k;
    int m, n, v, b;
    double s, s1, s2;
    double *dat = NULL;
    MAPTYPE *map = NULL;
    static double M[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    
    if (nrhs < 1 || nrhs > 3 || nlhs > 1)
    {
        mexErrMsgTxt("Incorrect usage.");
    }
    if (nrhs < 2)
    {
        s = mxGetNaN();
    }
    else
    {
        s = mxGetScalar(prhs[1]);
    }
    b = (nrhs < 3) ? 1 : mxGetScalar(prhs[2]);
    
    map = get_maps(prhs[0], &n);
    
    for(v=1; v<n; v++)
    {
        if (map[v].dim[0] != map[0].dim[0] ||
            map[v].dim[1] != map[0].dim[1] ||
            map[v].dim[2] != map[0].dim[2])
        {
            free_maps(map, n);
            mexErrMsgTxt("Incompatible image dimensions.");
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    
    k   = (map[0].dim[0])*(map[0].dim[1]);
    dat = (double *)mxCalloc(k, sizeof(double));
    
    for (v=0; v<n; v++)
    {
        if (mxIsNaN(s))
        {
            s1 = 0.0;
            m  = 0;
            for (i=0; i<map[v].dim[2]; i++)
            {
                M[14] = i+1;
                slice(M, dat, map[v].dim[0], map[v].dim[1], &map[v], 0, 0);
                for(j=0; j<k; j++)
                    if (mxIsFinite(dat[j]))
                    {
                        s1 += dat[j];
                        m++;
                    }
            }
            s1 /= (8.0*m);
        }
        else
        {
            s1 = s;
        }
        
        s2 = 0.0;
        m  = 0;
        for (i=0; i<map[v].dim[2]; i++)
        {
            M[14] = i+1;
            slice(M, dat, map[v].dim[0] ,map[v].dim[1], &map[v], 0, 0);
            for(j=0; j<k; j++)
                if (mxIsFinite(dat[j]) && dat[j]>s1)
                {
                    s2 += dat[j];
                    m++;
                }
        }
        if (b) s2 /= m;
        
        mxGetPr(plhs[0])[v] = s2;
    }
    
    free_maps(map, n);
    mxFree(dat);
}
