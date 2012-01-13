/*
 * $Id: spm_bsplinc.c 4624 2012-01-13 13:27:08Z john $
 * John Ashburner
 */
 
/*
 * This code is a modified version of that of Philippe Thevenaz, which I took from:
 *  http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs. Many thanks to Philippe Thevenaz for advice with the code.
 *
 * See:
 *  M. Unser, A. Aldroubi and M. Eden.
 *  "B-Spline Signal Processing: Part I-Theory,"
 *  IEEE Transactions on Signal Processing 41(2):821-832 (1993).
 *
 *  M. Unser, A. Aldroubi and M. Eden.
 *  "B-Spline Signal Processing: Part II-Efficient Design and Applications,"
 *  IEEE Transactions on Signal Processing 41(2):834-848 (1993).
 *
 *  M. Unser.
 *  "Splines: A Perfect Fit for Signal and Image Processing,"
 *  IEEE Signal Processing Magazine 16(6):22-38 (1999).
 *
*/

#include <math.h>
#include "mex.h"
#include "spm_mapping.h"
#include "bsplines.h"


/***************************************************************************************
Deconvolve the B-spline basis functions from the image volume
    vol - a handle for the volume to deconvolve
    c - the coefficients (arising from the deconvolution)
    d - the spline degree
    splinc0, splinc1, splinc2   - functions for 1D deconvolutions
*/
static int vol_coeffs(MAPTYPE *vol, double c[], int d[], void (*splinc[])())
{
    double  p[4];
    double *cp;
    int np;
    int i, j, k, n;
    double f[10240];

    /* Check that dimensions don't exceed size of f */
    if (vol->dim[1]>10240 ||vol->dim[2]>10240)
        return(1);

    /* Do a straight copy */
    cp = c;
    for(k=0; k<vol->dim[2]; k++)
    {
        double dk = k+1;
        for(j=0; j<vol->dim[1]; j++)
        {
            double dj = j+1;
            for(i=0;i<vol->dim[0];i++, cp++)
            {
                double di = i+1;
                resample(1,vol,cp,&di,&dj,&dk,0, 0.0);

                /* Not sure how best to handle NaNs */
                if (!mxIsFinite(*cp)) *cp = 0.0;
            }
        }
    }

    /* Deconvolve along the fastest dimension (X) */
    if (d[0]>1 && vol->dim[0]>1)
    {
        if (get_poles(d[0], &np, p)) return(1);
        for(k=0; k<vol->dim[2]; k++)
        {
            /* double dk = k+1; */
            for(j=0; j<vol->dim[1]; j++)
            {
                cp = &c[vol->dim[0]*(j+vol->dim[1]*k)];
                splinc[0](cp, vol->dim[0], p, np);
            }
        }
    }

    /* Deconvolve along the middle dimension (Y) */
    if (d[1]>1 && vol->dim[1]>1)
    {
        if (get_poles(d[1], &np, p)) return(1);
        n =vol->dim[0];
        for(k=0; k<vol->dim[2]; k++)
        {
            for(i=0;i<vol->dim[0];i++)
            {
                cp = &c[i+vol->dim[0]*vol->dim[1]*k];
                for(j=0; j<vol->dim[1]; j++, cp+=n)
                    f[j] = *cp;
                splinc[1](f, vol->dim[1], p, np);
                cp = &c[i+vol->dim[0]*vol->dim[1]*k];
                for(j=0; j<vol->dim[1]; j++, cp+=n)
                    *cp = f[j];
            }
        }
    }

    /* Deconvolve along the slowest dimension (Z) */
    if (d[2]>1 && vol->dim[2]>1)
    {
        if (get_poles(d[2], &np, p)) return(1);
        n = vol->dim[0]*vol->dim[1];
        for(j=0; j<vol->dim[1]; j++)
        {
            for(i=0;i<vol->dim[0];i++)
            {
                cp = &c[i+vol->dim[0]*j];
                for(k=0; k<vol->dim[2]; k++, cp+=n)
                    f[k] = *cp;
                splinc[2](f, vol->dim[2], p, np);
                cp = &c[i+vol->dim[0]*j];
                for(k=0; k<vol->dim[2]; k++, cp+=n)
                    *cp = f[k];
            }
        }
    }
    return(0);
}

/***************************************************************************************
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int k, d[3], sts;
    MAPTYPE *vol, *get_maps();
    double *c;
    void (*splinc[3])();

    if (nrhs < 2 || nlhs > 1)
        mexErrMsgTxt("Incorrect usage.");
    if (mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) ||
        (mxGetM(prhs[1])*mxGetN(prhs[1]) != 3 && mxGetM(prhs[1])*mxGetN(prhs[1]) != 6))
        mexErrMsgTxt("Incorrect usage.");

    for(k=0; k<3; k++)
    {
        d[k] = floor(mxGetPr(prhs[1])[k]+0.5);
        if (d[k]<0 || d[k]>7)
            mexErrMsgTxt("Bad spline degree.");
    }

    for(k=0; k<3; k++) splinc[k] = splinc_mirror;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) == 6)
    {
        for(k=0; k<3; k++)
            if (mxGetPr(prhs[1])[k+3])
                splinc[k] = splinc_wrap;
    }

    vol=get_maps(prhs[0], &k);
    if (k!=1)
    {
        free_maps(vol, k);
        mexErrMsgTxt("Too many images.");
    }

    plhs[0] = mxCreateNumericArray(3,vol->dim, mxDOUBLE_CLASS, mxREAL);
    c = (double *)mxGetPr(plhs[0]);

    sts = vol_coeffs(vol, c, d, splinc);

    if (sts)
    {
        free_maps(vol, k);
        mexErrMsgTxt("Problem with deconvolution.");
    }
    free_maps(vol, k);
}
