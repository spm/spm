/*
 * $Id: spm_render_vol.c 4453 2011-09-02 10:47:25Z guillaume $
 * John Ashburner
 */

#include <math.h>
#include "mex.h"
#include "spm_mapping.h"
#define RINT(A) floor((A)+0.5)

void surface(double mat[], double zbuff[], double xcords[], double ycords[], 
        double zcords[], mwSize xdim1, mwSize ydim1, MAPTYPE *vol, double thresh)
{
    /*
    Project the coordinates of voxels which lie on the surface of the object
    onto the viewing plane.
    A surface is defined as a voxel above "thresh" adjascent to voxel(s) below "thresh"
    - assuming 6 nearest neighbours.
    Only the voxels closest to the viewing plane are stored.

    Projection performed by matrix "mat".
    The elements of "mat" are numbered:
        0  4  8 12
        1  5  9 13
        2  6 10 14
        3  7 11 15

    voxel (x1, y1, z1) is projected onto viewing plane at:
    ((x1*mat[0] + y1*mat[4] + z1*mat[8] + mat[12]),
    (x1*mat[1] + y1*mat[5] + z1*mat[9] + mat[13]))
    Distance from the viewing plane (stored in "zbuff") is:
    (x1*mat[2] + y1*mat[6] + z1*mat[10] + mat[14])

    Note that elements 12-15 are unused, and that the projection is a parallel one.
    
    */
    int z, i,off1;
    double xs, ys;
    char *tmp, *current, *prev, *next;
    double mat1[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    double *dat;

    /* size of rectangle that is projected */
    xs = (0.01+sqrt(mat[0]*mat[0]+mat[4]*mat[4]+mat[8]*mat[8]));
    ys = (0.01+sqrt(mat[1]*mat[1]+mat[5]*mat[5]+mat[9]*mat[9]));

    dat     = (double *)mxCalloc(vol->dim[1]*vol->dim[0],sizeof(double));
    current = (char   *)mxCalloc(vol->dim[1]*vol->dim[0],sizeof(char));
    prev    = (char   *)mxCalloc(vol->dim[1]*vol->dim[0],sizeof(char));
    next    = (char   *)mxCalloc(vol->dim[1]*vol->dim[0],sizeof(char));

    mat1[14]=1.0;
    slice(mat1, dat, vol->dim[0],vol->dim[1], vol, 0,0);
    for(i=0; i<vol->dim[0]*vol->dim[1]; i++)
    {
        current[i] = ((dat[i]>thresh) ? 1 : 0);
        prev[i]=0;
    }

    for(z=1; z<=vol->dim[2]; z++)
    {
        int y;
        double x2, y2, z2;
        x2 = mat[12] + z*mat[8] -1.0;
        y2 = mat[13] + z*mat[9] -1.0;
        z2 = mat[14] + z*mat[10];

        mat1[14]=z+1.0;
        slice(mat1, dat, vol->dim[0],vol->dim[1], vol, 0,0);
        for(i=0; i<vol->dim[0]*vol->dim[1]; i++)
            next[i] = ((dat[i]>thresh) ? 1 : 0);
        off1 = 0;
        for(y=1; y<=vol->dim[1]; y++)
        {
            int x;
            double x3 = x2 + y*mat[4];
            double y3 = y2 + y*mat[5];
            double z3 = z2 + y*mat[6];

            for(x=1; x<=vol->dim[0]; x++)
            {
                if (current[off1] && (
                    x<=1 || x>=vol->dim[0] || !current[off1-          1] || !current[off1+          1] ||
                    y<=1 || y>=vol->dim[1] || !current[off1-vol->dim[0]] || !current[off1+vol->dim[0]] ||
                    !next[off1] || !prev[off1]))
                {
                    int ix4, iy4, ixstart, ixend, iystart, iyend;
                    double x4 = (x3 + x*mat[0]);
                    double y4 = (y3 + x*mat[1]);
                    double z4 = (z3 + x*mat[2]);
                    ixstart = (int)RINT(x4);
                    if (ixstart < 0) ixstart = 0;
                    ixend = (int)RINT(x4+xs);
                    if (ixend >= xdim1) ixend = xdim1-1;
                    iystart = (int)RINT(y4);
                    if (iystart < 0) iystart = 0;
                    iyend = (int)RINT(y4+ys);
                    if (iyend >= ydim1) iyend = ydim1-1;
                    for(iy4 = iystart ; iy4<=iyend; iy4++)
                        for(ix4 = ixstart ; ix4<=ixend; ix4++)
                        {
                            int off =  ix4+iy4*xdim1;
                            /* is new voxel closer than original one */
                            if (z4 < zbuff[off])
                            {
                                 zbuff[off] = z4;
                                xcords[off] = x;
                                ycords[off] = y;
                                zcords[off] = z;
                            }
                        }
                }
                off1++;
            }
        }
        tmp     = prev;
        prev    = current;
        current = next;
        next    = tmp;
    }
    mxFree(dat);
    mxFree(current);
    mxFree(prev);
    mxFree(next);
}

void render(double xcords[], double ycords[], double zcords[], mwSize xdim1, 
        mwSize ydim1, double out[], MAPTYPE *vol, double light[], double thresh, int nn)
{
    double *filt, *filtp;
    int dx, dy, dz, i, nnn;
    double *X, *Y, *Z, *dat;

    nnn  = (nn*2+1)*(nn*2+1)*(nn*2+1);
    X    = (double *)mxCalloc(nnn, sizeof(double));
    Y    = (double *)mxCalloc(nnn, sizeof(double));
    Z    = (double *)mxCalloc(nnn, sizeof(double));
    dat  = (double *)mxCalloc(nnn, sizeof(double));
    filt = (double *)mxCalloc((nn*2+1)*(nn*2+1)*(nn*2+1), sizeof(double));

    if (filt == NULL) return;

    /* create a 3D 'filter' which should produce the appropriate shading */
    filtp = filt;
    for(dz=-nn;dz<=nn;dz++)
            for(dy=-nn;dy<=nn;dy++)
                for(dx=-nn;dx<=nn;dx++)
                    *(filtp++) = 
                        (-light[0]*dx-light[1]*dy-light[2]*dz)/(nn*nn*nn*8);

    for(i=0; i<xdim1*ydim1; i++)
        if (xcords[i])
        {
            int x, y, z, j;
            double val;
            x = xcords[i];
            y = ycords[i];
            z = zcords[i];
            j = 0;
            for(dz=z-nn; dz<=z+nn; dz++)
                for(dy=y-nn; dy<=y+nn; dy++)
                    for(dx=x-nn; dx<=x+nn; dx++)
                    {
                        X[j] = dx;
                        Y[j] = dy;
                        Z[j] = dz;
                        j++;
                    }
            resample(nnn,vol,dat,X,Y,Z,0, 0.0);

            val = 0.0;
            for(j=0; j<nnn; j++)
            {
                if (dat[j]<thresh)
                    val+=filt[j];
            }
            if (val<0) val = 0.0;
            out[i] = val;
        }

    mxFree(filt);
    mxFree(X);
    mxFree(Y);
    mxFree(Z);
    mxFree(dat);
}

void initdat(double dat[], mwSize size, double val)
{
    mwIndex i;
    for(i=0; i<size; dat[i++] = val);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    MAPTYPE *map, *get_maps();
    mwSize m, n, k;
    int nn;
    double *mat;
    double *out, *odims, *params, *zbuff, *x, *y, *z;
    double light[3];

    if (nrhs != 4 || nlhs > 5) mexErrMsgTxt("Incorrect usage.");

    map = get_maps(prhs[0], &nn);
    if (nn!=1)
    {
        free_maps(map, nn);
        mexErrMsgTxt("Incorrect usage.");
    }

    for(k=1; k<nrhs; k++)
    {
        if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
            mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
        {
            free_maps(map, 1);
            mexErrMsgTxt("Arguments must be numeric, real, full and double.");
        }
    }

    /* get transformation matrix */
    if (mxGetM(prhs[1]) != 4 && mxGetN(prhs[1]) != 4)
    {
        free_maps(map, 1);
        mexErrMsgTxt("Transformation matrix must be 4 x 4.");
    }
    mat = mxGetPr(prhs[1]);

    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 2)
    {
        free_maps(map, 1);
        mexErrMsgTxt("Output dimensions must have two elements.");
    }
    odims = mxGetPr(prhs[2]);
    m = (mwSize)odims[0];
    n = (mwSize)odims[1];

    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 2)
    {
        free_maps(map, 1);
        mexErrMsgTxt("Parameters must have 2 elements.");
    }
    params = mxGetPr(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
    out = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(m,n, mxREAL);
    zbuff = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(m,n, mxREAL);
    x = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(m,n, mxREAL);
    y = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(m,n, mxREAL);
    z = mxGetPr(plhs[4]);

    initdat(zbuff,m*n,1024.0);
    surface(mat, zbuff, x, y, z, m, n, map, params[0]);

    light[0] = mat[2];
    light[1] = mat[6];
    light[2] = mat[10];

    render(x, y, z, m, n, out, map, light, params[0], (int)params[1]);
    free_maps(map, 1);
}
