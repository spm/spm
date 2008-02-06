/*
 * $Id: spm_affdef.c 1140 2008-02-06 19:24:05Z spm $
 * John Ashburner
 */

/* Note that according to the Matlab documentation, one should "avoid
   modifying input arguments in MEX-files".
   "In MATLAB 5.1 to 5.3.1, MATLAB arrays can share data.  There is
    currently no way for a MEX-file to determine that an array
    contains shared data.  MEX-files that modify their input arguments
    may corrupt arrays in the MATLAB workspace.  This style of programming
    is strongly discouraged."

   I have used this style of programming here in order to save memory.
*/
#include "mex.h"

static float *get_volume(const mxArray *ptr, int dims[3])
{
    int nd, i;
    const int *ldims;
    if (mxIsStruct(ptr) || !mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsSparse(ptr) || !mxIsSingle(ptr))
        mexErrMsgTxt("Data must be a single precision floating point multi-dimensional array.");

    nd = mxGetNumberOfDimensions(ptr);
    if (nd>3)
        mexErrMsgTxt("Too many dimensions in data.");

    ldims = mxGetDimensions(ptr);
    for(i=0; i<nd; i++)
        dims[i] = ldims[i];
    for(i=nd; i<3; i++)
        dims[i] = 1;

    return((float *)mxGetPr(ptr));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float  *y0=0,  *y1=0,  *y2=0, *p=0;
    int dim_g[3], dim_f[3];
    double *M;

        if (nrhs != 4 || nlhs >0)
                mexErrMsgTxt("Inappropriate usage.");

    y0 = get_volume(prhs[0], dim_g);
    y1 = get_volume(prhs[1], dim_f);
    if (dim_g[0] != dim_f[0] || dim_g[1] != dim_f[1] || dim_g[2] != dim_f[2])
        mexErrMsgTxt("Incompatible dimensions.");
    y2 = get_volume(prhs[2], dim_f);
    if (dim_g[0] != dim_f[0] || dim_g[1] != dim_f[1] || dim_g[2] != dim_f[2])
        mexErrMsgTxt("Incompatible dimensions.");

    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
        mxIsComplex(prhs[3]) || !mxIsDouble(prhs[3]) || mxGetM(prhs[3]) != 4 || mxGetN(prhs[3]) != 4)
        mexErrMsgTxt("Affine transform matrix must be 4x4.");

    M = mxGetPr(prhs[3]);

    p = y0+dim_g[0]*dim_g[1]*dim_g[2];

    while(y0<p)
    {
        float x0, x1, x2, x3;
        x0      = *y0;
        x1      = *y1;
        x2      = *y2;
        x3      =  M[3 + 0*4]*x0 + M[3 + 1*4]*x1 + M[3 + 2*4]*x2 + M[3 + 3*4];
        *(y0++) = (M[0 + 0*4]*x0 + M[0 + 1*4]*x1 + M[0 + 2*4]*x2 + M[0 + 3*4])/x3;
        *(y1++) = (M[1 + 0*4]*x0 + M[1 + 1*4]*x1 + M[1 + 2*4]*x2 + M[1 + 3*4])/x3;
        *(y2++) = (M[2 + 0*4]*x0 + M[2 + 1*4]*x1 + M[2 + 2*4]*x2 + M[2 + 3*4])/x3;
    }
}
