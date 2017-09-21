/* $Id: spm_mrf.c 7172 2017-09-21 16:31:30Z john $ */
/* (c) John Ashburner (2010) */

#include "mex.h"
#include <math.h>
#define MAXCLASSES 1024


static void mrf1(mwSize dm[], unsigned char q[], float p[], float G[], float w[], int code)
{
    mwSize i0, i1, i2, k, m, n;
    float a[MAXCLASSES], e[MAXCLASSES], *p0 = NULL, *p1 = NULL;
    unsigned char *q0 = NULL, *q1 = NULL;
    int it;

    m = dm[0]*dm[1]*dm[2];

    /* Use a red-black scheme, so the updates are for
       alternating voxels.  Then do another pass to
       update the other half.
       A B A B A B
       B A B A B A
       A B A B A B
       B A B A B A

       Updates involve computing the number of neighbours
       of each type (stored in vector a), and using the
       connectivity matrix (G) to update with:
           q = (p.*exp(G'*a))/sum((p.*exp(G'*a))
    */
    for(it=0; it<2; it++) 
    {
        mwSize i2start = it%2;
        for(i2=0; i2<dm[2]; i2++) /* Inferior -> Superior */
        {
            mwSize i1start = (i2start == (i2%2));
            for(i1=0; i1<dm[1]; i1++) /* Posterior -> Anterior */
            {
                mwSize i0start = (i1start == (i1%2));
                p1 = p + dm[0]*(i1+dm[1]*i2);
                q1 = q + dm[0]*(i1+dm[1]*i2);

                for(i0=i0start; i0<dm[0]; i0+=2) /* Left -> Right */
                {
                    float se, mx;
                    unsigned char *qq = NULL;

                    /* Pointers to current voxel in first volume */
                    p0 = p1 + i0;
                    q0 = q1 + i0;

                    /* Initialise neighbour counts to zero */
                    for(k=0; k<dm[3]; k++) a[k] = 0.0;

                    /* Count neighbours of each class */
                    if(i2>0)       /* Inferior */
                    {

                        qq = q0 - dm[0]*dm[1];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[2];
                    }

                    if(i2<dm[2]-1) /* Superior */
                    {
                        qq = q0 + dm[0]*dm[1];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[2];
                    }

                    if(i1>0)       /* Posterior */
                    {
                        qq = q0 - dm[0];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[1];
                    }

                    if(i1<dm[1]-1) /* Anterior */
                    {
                        qq = q0 + dm[0];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[1];
                    }

                    if(i0>0)       /* Left */
                    {
                        qq = q0 - 1;
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[0];
                    }

                    if(i0<dm[0]-1) /* Right */
                    {
                        qq = q0 + 1;
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[0];
                    }

                    /* Responsibility data is uint8, so correct scaling.
                       Note also that data is divided by 6 (the number
                       of neighbours examined). */
                    for(k=0; k<dm[3]; k++)
                        a[k]/=(255.0*6.0);

                    if (code == 1) 
                    {
                        /* Weights are in the form of a matrix,
                           shared among all voxels. */
                        float *g;
                        for(k=0, g=G; k<dm[3]; k++)
                        {
                            e[k] = 0;
                            for(n=0; n<dm[3]; n++, g++)
                                e[k] += (*g)*a[n];
                        }
                    }
                    else if (code == 2)
                    {
                        /* Weights are assumed to be a diagonal matrix,
                           so only the diagonal elements are passed. */
                        for(k=0; k<dm[3]; k++)
                            e[k] = G[k]*a[k];
                    }
                    else if (code == 3)
                    {
                        /* Separate weights for each voxel, in the form of
                           the full matrix (loads of memory). */
                        float *g;
                        g = G + i0+dm[0]*(i1+dm[1]*i2);
                        for(k=0; k<dm[3]; k++)
                        {
                            e[k] = 0.0;
                            for(n=0; n<dm[3]; n++, g+=m)
                                e[k] += (*g)*a[n];
                        }
                    }
                    else if (code == 4)
                    {
                        /* Separate weight matrices for each voxel,
                           where the matrices are assumed to be symmetric
                           with zeros on the diagonal. For a 4x4
                           matrix, the elements are ordered as
                           (2,1), (3,1), (4,1), (3,2), (4,2), (4,3).
                         */
                        float *g;
                        g = G + i0+dm[0]*(i1+dm[1]*i2);
                        for(k=0; k<dm[3]; k++) e[k] = 0.0;

                        for(k=0; k<dm[3]; k++)
                        {
                            for(n=k+1; n<dm[3]; n++, g+=m)
                            {
                                e[k] += (*g)*a[n];
                                e[n] += (*g)*a[k];
                            }
                        }
                    }
                    else if (code == 5)
                    {
                        /* Separate weight matrices for each voxel,
                           where the matrices are assumed to be symmetric
                           with zeros on the diagonal. For a 4x4
                           matrix, the elements are ordered as
                           (2,1), (3,1), (4,1), (3,2), (4,2), (4,3).
                           
                           The weight matrices are encoded as uint8, and
                           their values need to be scaled by -0.0625 to
                           bring them into a reasonable range.
                         */
                        unsigned char *g;
                        g = (unsigned char *)G + i0+dm[0]*(i1+dm[1]*i2);
                        for(k=0; k<dm[3]; k++) e[k] = 0.0;

                        for(k=0; k<dm[3]; k++)
                        {
                            for(n=k+1; n<dm[3]; n++, g+=m)
                            {
                                e[k] += ((float)(*g))*a[n];
                                e[n] += ((float)(*g))*a[k];
                            }
                        }
                        for(k=0; k<dm[3]; k++)
                            e[k] = -0.0625*e[k];
                    }

                    /* Prevent overflow by subtracting maximum value before computing exp - based on log-sum-exp trick */
                    mx = -3.4028e+38;
                    for(k=0; k<dm[3]; k++)
                        if (e[k]>mx) mx = e[k];

                    se = 0.0;
                    for(k=0; k<dm[3]; k++)
                    {
                        e[k] = exp((double)(e[k]-mx))*p0[k*m];
                        se  += e[k];
                    }

                    /* Normalise responsibilities to sum to 1
                       and rescale for saving as uint8 data. */
                    se = 255.0/se;
                    for(k=0; k<dm[3]; k++)
                        q0[k*m] = (unsigned char)(e[k]*se+0.5);
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize i;
    mwSize dm[4];
    float *p = NULL, w[3];
    unsigned char *q = NULL;
    float *G = NULL;
    int code=0;

    if (nrhs<3 || nrhs>4 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsUint8(prhs[0]))
        mexErrMsgTxt("First arg must be numeric, real, full and uint8.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsSingle(prhs[1]))
        mexErrMsgTxt("Second arg must be numeric, real, full and single.");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]))
        mexErrMsgTxt("Third arg must be numeric, real and full.");

    if (mxGetNumberOfDimensions(prhs[0])!= mxGetNumberOfDimensions(prhs[1]) ||
        mxGetNumberOfDimensions(prhs[0])>4)
        mexErrMsgTxt("First or second args have wrong number of dimensions.");

    for(i=0; i<mxGetNumberOfDimensions(prhs[0]); i++)
        dm[i] = mxGetDimensions(prhs[0])[i];

    for(i=mxGetNumberOfDimensions(prhs[0]); i<4; i++)
        dm[i] = 1;

    if (dm[3]>MAXCLASSES) mexErrMsgTxt("Too many classes.");

    for(i=0; i<4; i++)
        if (mxGetDimensions(prhs[1])[i] != dm[i])
            mexErrMsgTxt("First and second args have incompatible dimensions.");

    if (mxGetDimensions(prhs[2])[1] == 1)
    {
        code = 2;
        if (mxGetDimensions(prhs[2])[0] != dm[3])
            mexErrMsgTxt("Third arg has incompatible dimensions.");

        if (!mxIsSingle(prhs[2])) mexErrMsgTxt("Third arg must be single.");
    }
    else if (mxGetNumberOfDimensions(prhs[2])==2)
    {
        code = 1;
        if (mxGetDimensions(prhs[2])[0] != dm[3] || mxGetDimensions(prhs[2])[1] != dm[3])
            mexErrMsgTxt("Third arg has incompatible dimensions.");

        if (!mxIsSingle(prhs[2])) mexErrMsgTxt("Third arg must be single.");
    }
    else if (mxGetNumberOfDimensions(prhs[2])==5)
    {
        code = 3;
        for(i=0; i<4; i++)
            if (mxGetDimensions(prhs[2])[i] != dm[i])
                mexErrMsgTxt("Third arg has incompatible dimensions.");

        if (mxGetDimensions(prhs[2])[4] != dm[3])
            mexErrMsgTxt("Third arg has incompatible dimensions.");

        if (!mxIsSingle(prhs[2])) mexErrMsgTxt("Third arg must be single.");
    }
    else if (mxGetNumberOfDimensions(prhs[2])==4)
    {
        for(i=0; i<3; i++)
            if (mxGetDimensions(prhs[2])[i] != dm[i])
                mexErrMsgTxt("Third arg has incompatible dimensions.");

        if (mxGetDimensions(prhs[2])[3] != (dm[3]*(dm[3]-1))/2)
            mexErrMsgTxt("Third arg has incompatible dimensions.");

        if (mxIsSingle(prhs[2]))
            code = 4;
        else if (mxIsUint8(prhs[0]))
            code = 5;
        else
            mexErrMsgTxt("Third arg must be either single or uint8.");
    }
    else
        mexErrMsgTxt("Third arg has incompatible dimensions.");


    p = (float *)mxGetData(prhs[1]);
    G = (float *)mxGetData(prhs[2]);

    if (nrhs>=4)
    {
        /* Adjustment for anisotropic voxel sizes.  w should contain
           the square of each voxel size. */
        if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsSingle(prhs[3]))
            mexErrMsgTxt("Fourth arg must be numeric, real, full and single.");

        if (mxGetNumberOfElements(prhs[3]) != 3)
            mexErrMsgTxt("Fourth arg must contain three elements.");

        for(i=0; i<3; i++) w[i] = ((float *)mxGetData(prhs[3]))[i];
    }
    else
    {
        for(i=0; i<3; i++) w[i] = 1.0;
    }

    if (nlhs>0)
    {
        /* Copy input to output */
        unsigned char *q0;
        plhs[0]  = mxCreateNumericArray(4,dm, mxUINT8_CLASS, mxREAL);
        q0 = (unsigned char *)mxGetData(prhs[0]);
        q  = (unsigned char *)mxGetData(plhs[0]);

        for(i=0; i<dm[0]*dm[1]*dm[2]*dm[3]; i++)
            q[i] = q0[i];
    }
    else /* Note the nasty side effects - but it does save memory */
        q = (unsigned char *)mxGetData(prhs[0]);

    mrf1(dm, q,p,G,w,code);
}

