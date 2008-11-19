/*
 * $Id: spm_warp.c 2480 2008-11-19 17:47:49Z john $
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

#include <math.h>
#include <stdio.h>
#include "mex.h"

static float resample(unsigned char vol[], float x1, float x2, float x3, int dim[3])
{
    float out;
    if (x3>=1 && x3<dim[2] && x2>=1 && x2<dim[1] && x1>=1 && x1<dim[0])
    {
        float k111,k112,k121,k122,k211,k212,k221,k222;
        float dx11, dx12, dx21, dx22, dx31, dx32;
        int off1, off2, x1coord, x2coord, x3coord;

        x1coord = (int)floor(x1); dx11=x1-x1coord; dx12=1.0-dx11;
        x2coord = (int)floor(x2); dx21=x2-x2coord; dx22=1.0-dx21;
        x3coord = (int)floor(x3); dx31=x3-x3coord; dx32=1.0-dx31;
 
        off1 = x1coord-1 + dim[0]*(x2coord-1 + dim[1]*(x3coord-1));
        k222 = vol[off1]; k122 = vol[off1+1]; off2 = off1+dim[0];
        k212 = vol[off2]; k112 = vol[off2+1]; off1+= dim[0]*dim[1];
        k221 = vol[off1]; k121 = vol[off1+1]; off2 = off1+dim[0];
        k211 = vol[off2]; k111 = vol[off2+1];

        out =  (((k222*dx12 + k122*dx11)*dx22 + (k212*dx12 + k112*dx11)*dx21))*dx32
             + (((k221*dx12 + k121*dx11)*dx22 + (k211*dx12 + k111*dx11)*dx21))*dx31;
    }
    else out = mxGetNaN();
    return(out);
}

static float resample_d(unsigned char vol[], float x1, float x2, float x3, int dim[3], float *gradx1, float *gradx2, float *gradx3)
{
    float out;
    if (x3>=1 && x3<dim[2] && x2>=1 && x2<dim[1] && x1>=1 && x1<dim[0])
    {
        float k111,k112,k121,k122,k211,k212,k221,k222;
        float dx11, dx12, dx21, dx22, dx31, dx32;
        int off1, off2, x1coord, x2coord, x3coord;

        x1coord = (int)floor(x1); dx11=x1-x1coord; dx12=1.0-dx11;
        x2coord = (int)floor(x2); dx21=x2-x2coord; dx22=1.0-dx21;
        x3coord = (int)floor(x3); dx31=x3-x3coord; dx32=1.0-dx31;
 
        off1 = x1coord-1  + dim[0]*(x2coord-1  + dim[1]*(x3coord-1));
        k222 = vol[off1]; k122 = vol[off1+1]; off2 = off1+dim[0];
        k212 = vol[off2]; k112 = vol[off2+1]; off1+= dim[0]*dim[1];
        k221 = vol[off1]; k121 = vol[off1+1]; off2 = off1+dim[0];
        k211 = vol[off2]; k111 = vol[off2+1];

        *gradx1 = (((k111 - k211)*dx21 + (k121 - k221)*dx22))*dx31
                + (((k112 - k212)*dx21 + (k122 - k222)*dx22))*dx32;

        k111 = (k111*dx11 + k211*dx12);
        k121 = (k121*dx11 + k221*dx12);
        k112 = (k112*dx11 + k212*dx12);
        k122 = (k122*dx11 + k222*dx12);

        *gradx2 = (k111 - k121)*dx31 + (k112 - k122)*dx32;

        k111 = k111*dx21 + k121*dx22;
        k112 = k112*dx21 + k122*dx22;

        *gradx3 = k111 - k112;
        out     = k111*dx31 + k112*dx32;
    }
    else
    {
        out     = mxGetNaN();
        *gradx1 = 0.0;
        *gradx2 = 0.0;
        *gradx3 = 0.0;
    }
    return(out);
}

#ifndef NOZEROMASK
#define known(x) (((x)!=0) && mxIsFinite(x))
#endif
#ifdef NOZEROMASK
#define known(x) mxIsFinite(x)
#endif

static void invertX(float X[4][3], float IX[4][4], float *dt)
{
    float id;
    *dt = X[0][0]*(X[3][1]*(X[1][2]-X[2][2])+X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2]))+
          X[1][0]*(X[3][2]*(X[0][1]-X[2][1])+X[0][2]*(X[2][1]-X[3][1])+X[2][2]*(X[3][1]-X[0][1]))+
          X[2][0]*(X[0][1]*(X[1][2]-X[3][2])+X[3][1]*(X[0][2]-X[1][2])+X[1][1]*(X[3][2]-X[0][2]))+
          X[3][0]*(X[1][2]*(X[2][1]-X[0][1])+X[0][2]*(X[1][1]-X[2][1])+X[2][2]*(X[0][1]-X[1][1]));
    id = 1/(*dt);
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

static void getJ(float Y[4][3], float IX[4][4], float J[3][3])
{
    J[0][0] = Y[0][0]*IX[0][0] + Y[1][0]*IX[0][1] + Y[2][0]*IX[0][2] + Y[3][0]*IX[0][3];
    J[0][1] = Y[0][1]*IX[0][0] + Y[1][1]*IX[0][1] + Y[2][1]*IX[0][2] + Y[3][1]*IX[0][3];
    J[0][2] = Y[0][2]*IX[0][0] + Y[1][2]*IX[0][1] + Y[2][2]*IX[0][2] + Y[3][2]*IX[0][3];

    J[1][0] = Y[0][0]*IX[1][0] + Y[1][0]*IX[1][1] + Y[2][0]*IX[1][2] + Y[3][0]*IX[1][3];
    J[1][1] = Y[0][1]*IX[1][0] + Y[1][1]*IX[1][1] + Y[2][1]*IX[1][2] + Y[3][1]*IX[1][3];
    J[1][2] = Y[0][2]*IX[1][0] + Y[1][2]*IX[1][1] + Y[2][2]*IX[1][2] + Y[3][2]*IX[1][3];

    J[2][0] = Y[0][0]*IX[2][0] + Y[1][0]*IX[2][1] + Y[2][0]*IX[2][2] + Y[3][0]*IX[2][3];
    J[2][1] = Y[0][1]*IX[2][0] + Y[1][1]*IX[2][1] + Y[2][1]*IX[2][2] + Y[3][1]*IX[2][3];
    J[2][2] = Y[0][2]*IX[2][0] + Y[1][2]*IX[2][1] + Y[2][2]*IX[2][2] + Y[3][2]*IX[2][3];
}

static void invertJ(float J[3][3], float IJ[3][3], float *dt)
{
    float id;
    float t0, t1, t2;

    t0 = J[1][1]*J[2][2]-J[1][2]*J[2][1];
    t1 = J[1][2]*J[2][0]-J[1][0]*J[2][2];
    t2 = J[1][0]*J[2][1]-J[1][1]*J[2][0];
    *dt = J[0][0]*t0 + J[0][1]*t1 + J[0][2]*t2;
    id = 1.0/(*dt);

    IJ[0][0] = t0*id;
    IJ[0][1] = (J[0][2]*J[2][1]-J[0][1]*J[2][2])*id;
    IJ[0][2] = (J[0][1]*J[1][2]-J[0][2]*J[1][1])*id;

    IJ[1][0] = t1*id;
    IJ[1][1] = (J[0][0]*J[2][2]-J[0][2]*J[2][0])*id;
    IJ[1][2] = (J[0][2]*J[1][0]-J[0][0]*J[1][2])*id;

    IJ[2][0] = t2*id;
    IJ[2][1] = (J[0][1]*J[2][0]-J[0][0]*J[2][1])*id;
    IJ[2][2] = (J[0][0]*J[1][1]-J[0][1]*J[1][0])*id;
}

static void getdIJ(float J[3][3], float ix00, float ix10, float ix20,
    float IJ[3][3], float dIJ[3][3][3], float *dt, float *dt0, float *dt1, float *dt2)
{
    float id, id2, ddt;
    float t00,t10,t20,t01,t11,t21,t02,t12,t22;

    t00 = J[1][1]*J[2][2]-J[1][2]*J[2][1];
    t01 = J[0][2]*J[2][1]-J[0][1]*J[2][2];
    t02 = J[0][1]*J[1][2]-J[0][2]*J[1][1];

    t10 = J[1][2]*J[2][0]-J[1][0]*J[2][2];
    t11 = J[0][0]*J[2][2]-J[0][2]*J[2][0];
    t12 = J[0][2]*J[1][0]-J[0][0]*J[1][2];

    t20 = J[1][0]*J[2][1]-J[1][1]*J[2][0];
    t21 = J[0][1]*J[2][0]-J[0][0]*J[2][1];
    t22 = J[0][0]*J[1][1]-J[0][1]*J[1][0];

    *dt = J[0][0]*t00 + J[0][1]*t10 + J[0][2]*t20;
    id = 1.0/(*dt);

    IJ[0][0] = t00*id;  IJ[0][1] = t01*id;  IJ[0][2] = t02*id;
    IJ[1][0] = t10*id;  IJ[1][1] = t11*id;  IJ[1][2] = t12*id;
    IJ[2][0] = t20*id;  IJ[2][1] = t21*id;  IJ[2][2] = t22*id;

    id2 = id*id;

    *dt0 = ix00*t00 + ix10*t01 + ix20*t02;
    ddt = id2*(*dt0);
    dIJ[0][0][0] = -t00*ddt;
    dIJ[0][0][1] = -t01*ddt;
    dIJ[0][0][2] = -t02*ddt;
    dIJ[0][1][0] = (J[1][2]*ix20-ix10*J[2][2])*id - t10*ddt;
    dIJ[0][1][1] = (ix00*J[2][2]-J[0][2]*ix20)*id - t11*ddt;
    dIJ[0][1][2] = (J[0][2]*ix10-ix00*J[1][2])*id - t12*ddt;
    dIJ[0][2][0] = (ix10*J[2][1]-J[1][1]*ix20)*id - t20*ddt;
    dIJ[0][2][1] = (J[0][1]*ix20-ix00*J[2][1])*id - t21*ddt;
    dIJ[0][2][2] = (ix00*J[1][1]-J[0][1]*ix10)*id - t22*ddt;

    *dt1 = ix00*t10 + ix10*t11 + ix20*t12;
    ddt = id2*(*dt1);
    dIJ[1][0][0] = (ix10*J[2][2]-J[1][2]*ix20)*id - t00*ddt;
    dIJ[1][0][1] = (J[0][2]*ix20-ix00*J[2][2])*id - t01*ddt;
    dIJ[1][0][2] = (ix00*J[1][2]-J[0][2]*ix10)*id - t02*ddt;
    dIJ[1][1][0] = -t10*ddt;
    dIJ[1][1][1] = -t11*ddt;
    dIJ[1][1][2] = -t12*ddt;
    dIJ[1][2][0] = (J[1][0]*ix20-ix10*J[2][0])*id - t20*ddt;
    dIJ[1][2][1] = (ix00*J[2][0]-J[0][0]*ix20)*id - t21*ddt;
    dIJ[1][2][2] = (J[0][0]*ix10-ix00*J[1][0])*id - t22*ddt;

    *dt2 = ix00*t20 + ix10*t21 + ix20*t22;
    ddt = id2*(*dt2);
    dIJ[2][0][0] = (J[1][1]*ix20-ix10*J[2][1])*id - t00*ddt;
    dIJ[2][0][1] = (ix00*J[2][1]-J[0][1]*ix20)*id - t01*ddt;
    dIJ[2][0][2] = (J[0][1]*ix10-ix00*J[1][1])*id - t02*ddt;
    dIJ[2][1][0] = (ix10*J[2][0]-J[1][0]*ix20)*id - t10*ddt;
    dIJ[2][1][1] = (J[0][0]*ix20-ix00*J[2][0])*id - t11*ddt;
    dIJ[2][1][2] = (ix00*J[1][0]-J[0][0]*ix10)*id - t12*ddt;
    dIJ[2][2][0] = -t20*ddt;
    dIJ[2][2][1] = -t21*ddt;
    dIJ[2][2][2] = -t22*ddt;
}

static float sumJ2(float J[3][3])
{
    return(J[0][0]*J[0][0]+J[0][1]*J[0][1]+J[0][2]*J[0][2]+
           J[1][0]*J[1][0]+J[1][1]*J[1][1]+J[1][2]*J[1][2]+
           J[2][0]*J[2][0]+J[2][1]*J[2][1]+J[2][2]*J[2][2]);
}

static void derivsumJ2(float J[3][3], float ix00, float ix10, float ix20, float *d0, float *d1, float *d2)
{
    *d0 = 2*(J[0][0]*ix00 + J[1][0]*ix10 + J[2][0]*ix20);
    *d1 = 2*(J[0][1]*ix00 + J[1][1]*ix10 + J[2][1]*ix20);
    *d2 = 2*(J[0][2]*ix00 + J[1][2]*ix10 + J[2][2]*ix20);
}

static void derivsumIJ2(float IJ[3][3], float dIJ[3][3][3], float *d0, float *d1, float *d2)
{
    *d0 = 2*(IJ[0][0]*dIJ[0][0][0] + IJ[0][1]*dIJ[0][0][1] + IJ[0][2]*dIJ[0][0][2] +
         IJ[1][0]*dIJ[0][1][0] + IJ[1][1]*dIJ[0][1][1] + IJ[1][2]*dIJ[0][1][2] +
         IJ[2][0]*dIJ[0][2][0] + IJ[2][1]*dIJ[0][2][1] + IJ[2][2]*dIJ[0][2][2]);
    *d1 = 2*(IJ[0][0]*dIJ[1][0][0] + IJ[0][1]*dIJ[1][0][1] + IJ[0][2]*dIJ[1][0][2] +
         IJ[1][0]*dIJ[1][1][0] + IJ[1][1]*dIJ[1][1][1] + IJ[1][2]*dIJ[1][1][2] +
         IJ[2][0]*dIJ[1][2][0] + IJ[2][1]*dIJ[1][2][1] + IJ[2][2]*dIJ[1][2][2]);
    *d2 = 2*(IJ[0][0]*dIJ[2][0][0] + IJ[0][1]*dIJ[2][0][1] + IJ[0][2]*dIJ[2][0][2] +
         IJ[1][0]*dIJ[2][1][0] + IJ[1][1]*dIJ[2][1][1] + IJ[1][2]*dIJ[2][1][2] +
         IJ[2][0]*dIJ[2][2][0] + IJ[2][1]*dIJ[2][2][1] + IJ[2][2]*dIJ[2][2][2]);
}


typedef struct
{
    float ix[4][4];
    float volx;
    int   o1, o2, o3;
    unsigned long mask;
} Tettype;

static void fun_ss(Tettype *tet, float Y[4][3], float *hp)
{
    static float J[3][3], IJ[3][3];
    float dt;
    getJ(Y, tet->ix, J);
    invertJ(J, IJ, &dt);
    *hp += 0.25*tet->volx*(1+dt)*(sumJ2(J) + sumJ2(IJ) - 6.0);
}

static void dfun_ss(Tettype *tet, float Y[4][3],
          float *hp, float *dhp0, float *dhp1, float *dhp2)
{
    static float J[3][3], IJ[3][3], dIJ[3][3][3];
    float ss, dss0f, dss1f, dss2f, dss0b, dss1b, dss2b, dt, dd0, dd1, dd2;

    getJ(Y, tet->ix, J);
    getdIJ(J, tet->ix[0][0], tet->ix[1][0], tet->ix[2][0], IJ, dIJ, &dt, &dd0, &dd1, &dd2);
    ss = sumJ2(J) + sumJ2(IJ) - 6.0;
    derivsumJ2(J, tet->ix[0][0], tet->ix[1][0], tet->ix[2][0], &dss0f, &dss1f, &dss2f);
    derivsumIJ2(IJ, dIJ, &dss0b, &dss1b, &dss2b);

    *hp   += 0.25*tet->volx * (1+dt)*ss;
    *dhp0 += 0.25*tet->volx *((1+dt)*(dss0f+dss0b)+dd0*ss);
    *dhp1 += 0.25*tet->volx *((1+dt)*(dss1f+dss1b)+dd1*ss);
    *dhp2 += 0.25*tet->volx *((1+dt)*(dss2f+dss2b)+dd2*ss);
}

static Tettype oddtets[32], eventets[8];

static void setup_consts_sub(int dim[3], float x[3][3], Tettype *tet, float vox_g[3])
{
    int i;
    float X[4][3], dtx;
    X[0][0] = 0.0;
    X[0][1] = 0.0;
    X[0][2] = 0.0;

    X[1][0] = x[0][0]*vox_g[0];
    X[1][1] = x[0][1]*vox_g[1];
    X[1][2] = x[0][2]*vox_g[2];

    X[2][0] = x[1][0]*vox_g[0];
    X[2][1] = x[1][1]*vox_g[1];
    X[2][2] = x[1][2]*vox_g[2];

    X[3][0] = x[2][0]*vox_g[0];
    X[3][1] = x[2][1]*vox_g[1];
    X[3][2] = x[2][2]*vox_g[2];
    invertX(X, tet->ix, &dtx);
    tet->volx = fabs((double)dtx)/6.0;

    tet->o1 = x[0][0]+dim[0]*(x[0][1]+dim[1]*x[0][2]);
    tet->o2 = x[1][0]+dim[0]*(x[1][1]+dim[1]*x[1][2]);
    tet->o3 = x[2][0]+dim[0]*(x[2][1]+dim[1]*x[2][2]);

    tet->mask = 0;
    for(i=0; i<3; i++)
    {
        if(x[i][0] < 0) tet->mask |=  1;
        if(x[i][0] > 0) tet->mask |=  2;
        if(x[i][1] < 0) tet->mask |=  4;
        if(x[i][1] > 0) tet->mask |=  8;
        if(x[i][2] < 0) tet->mask |= 16;
        if(x[i][2] > 0) tet->mask |= 32;
    }
}

static void setup_consts(int dim[3], float vox_g[3])
{
    int i;
    static float xo[][3][3] = {
        {{ 1, 0, 0}, { 1, 0, 1}, { 1, 1, 0}},
        {{ 0, 1, 1}, { 1, 0, 1}, { 0, 0, 1}},
        {{ 0, 1, 1}, { 0, 1, 0}, { 1, 1, 0}},
        {{ 1, 1, 0}, { 1, 0, 1}, { 0, 1, 1}},

        {{-1, 0, 0}, {-1, 0, 1}, {-1,-1, 0}},
        {{ 0,-1, 1}, {-1, 0, 1}, { 0, 0, 1}},
        {{ 0,-1, 1}, { 0,-1, 0}, {-1,-1, 0}},
        {{-1,-1, 0}, {-1, 0, 1}, { 0,-1, 1}},

        {{-1, 0, 0}, {-1, 0,-1}, {-1, 1, 0}},
        {{ 0, 1,-1}, {-1, 0,-1}, { 0, 0,-1}},
        {{ 0, 1,-1}, { 0, 1, 0}, {-1, 1, 0}},
        {{-1, 1, 0}, {-1, 0,-1}, { 0, 1,-1}},

        {{ 1, 0, 0}, { 1, 0,-1}, { 1,-1, 0}},
        {{ 0,-1,-1}, { 1, 0,-1}, { 0, 0,-1}},
        {{ 0,-1,-1}, { 0,-1, 0}, { 1,-1, 0}},
        {{ 1,-1, 0}, { 1, 0,-1}, { 0,-1,-1}},


        {{-1, 0, 1}, {-1, 0, 0}, {-1, 1, 0}},
        {{-1, 0, 1}, { 0, 1, 1}, { 0, 0, 1}},
        {{ 0, 1, 0}, { 0, 1, 1}, {-1, 1, 0}},
        {{-1, 0, 1}, {-1, 1, 0}, { 0, 1, 1}},

        {{ 1, 0, 1}, { 1, 0, 0}, { 1,-1, 0}},
        {{ 1, 0, 1}, { 0,-1, 1}, { 0, 0, 1}},
        {{ 0,-1, 0}, { 0,-1, 1}, { 1,-1, 0}},
        {{ 1, 0, 1}, { 1,-1, 0}, { 0,-1, 1}},

        {{ 1, 0,-1}, { 1, 0, 0}, { 1, 1, 0}},
        {{ 1, 0,-1}, { 0, 1,-1}, { 0, 0,-1}},
        {{ 0, 1, 0}, { 0, 1,-1}, { 1, 1, 0}},
        {{ 1, 0,-1}, { 1, 1, 0}, { 0, 1,-1}},

        {{-1, 0,-1}, {-1, 0, 0}, {-1,-1, 0}},
        {{-1, 0,-1}, { 0,-1,-1}, { 0, 0,-1}},
        {{ 0,-1, 0}, { 0,-1,-1}, {-1,-1, 0}},
        {{-1, 0,-1}, {-1,-1, 0}, { 0,-1,-1}}};

    static float xe[][3][3] = {
        {{ 1, 0, 0}, { 0, 0, 1}, { 0, 1, 0}},
        {{-1, 0, 0}, { 0, 0, 1}, { 0,-1, 0}},
        {{-1, 0, 0}, { 0, 0,-1}, { 0, 1, 0}},
        {{ 1, 0, 0}, { 0, 0,-1}, { 0,-1, 0}},

        {{ 0, 0, 1}, {-1, 0, 0}, { 0, 1, 0}},
        {{ 0, 0, 1}, { 1, 0, 0}, { 0,-1, 0}},
        {{ 0, 0,-1}, { 1, 0, 0}, { 0, 1, 0}},
        {{ 0, 0,-1}, {-1, 0, 0}, { 0,-1, 0}}};

    for(i=0; i<32; i++) setup_consts_sub(dim, xo[i],  &oddtets[i], vox_g);
    for(i=0; i< 8; i++) setup_consts_sub(dim, xe[i], &eventets[i], vox_g);
}


static void tweek(int x0, int x1, int x2, float *Y0, float *Y1, float *Y2,
    int dim_g[3], unsigned char gvol[],
    int dim_f[3], float vox_f[3], unsigned char fvol[],
    int msk, float lambda, float epsilon, float v, float scale,
    float *hf, float *hp, int *n, float *sumf, float *sumg, int *cnt, int sgn)
{
    float h0, dh0, dh1, dh2, t0, t1, t2, tt0, tt1, tt2, d;
    float hf0, hf1;
    float hp0, hp1, dhp0, dhp1, dhp2;
    float f, df0, df1, df2, g;
    float *y0, *y1, *y2;
    int i, iter, o;
    int flg = 0;
    Tettype *tets;
    int ntets;

    flg = 0;
    if ((x0+x1+x2)%2){ tets =  oddtets; ntets = 32;}
    else             { tets = eventets; ntets = 8;}

    o   = x0-1 + dim_g[0]*(x1-1 + dim_g[1]*(x2-1));
    y0  = Y0+o; y1 = Y1+o; y2 = Y2+o;
    t0  = y0[0]; t1 = y1[0]; t2 = y2[0];
    hp0 = dhp0 = dhp1 = dhp2 = 0.0;

    for(i=0; i<ntets; i++)
    {
        if (!(tets[i].mask & msk))
        {
            float Y[4][3];
            Y[0][0] = t0*vox_f[0];             Y[0][1] = t1*vox_f[1];             Y[0][2] = t2*vox_f[2];
            Y[1][0] = y0[tets[i].o1]*vox_f[0]; Y[1][1] = y1[tets[i].o1]*vox_f[1]; Y[1][2] = y2[tets[i].o1]*vox_f[2];
            Y[2][0] = y0[tets[i].o2]*vox_f[0]; Y[2][1] = y1[tets[i].o2]*vox_f[1]; Y[2][2] = y2[tets[i].o2]*vox_f[2];
            Y[3][0] = y0[tets[i].o3]*vox_f[0]; Y[3][1] = y1[tets[i].o3]*vox_f[1]; Y[3][2] = y2[tets[i].o3]*vox_f[2];
            dfun_ss(&(tets[i]), Y, &hp0, &dhp0, &dhp1, &dhp2);
        }
    }
    dhp0 *= vox_f[0];
    dhp1 *= vox_f[1];
    dhp2 *= vox_f[2];

    f    = resample_d(fvol, t0, t1, t2, dim_f, &df0,&df1,&df2);
    if (!known(f))
    {
        flg  = 1;
        h0   =  hp0*lambda;
        dh0  = dhp0*lambda;
        dh1  = dhp1*lambda;
        dh2  = dhp2*lambda;
    }
    else
    {
        g      = gvol[o];
        d      = (f*scale-g);
        hf0    = d*d/2;
        d     *= scale/v;

        h0     = hf0/v+ hp0*lambda;
        dh0    = d*df0+dhp0*lambda;
        dh1    = d*df1+dhp1*lambda;
        dh2    = d*df2+dhp2*lambda;
    }

    if (!known(h0))
    {
        (void)mexPrintf("X");
        return;
    }

    /* Bracket such that determinant is positive */
    tt0 = t0 - dh0*epsilon;
    tt1 = t1 - dh1*epsilon;
    tt2 = t2 - dh2*epsilon;

    for(i=0; i<ntets; i++)
    {
        if (!(tets[i].mask & msk))
        {
            iter=0;
            while(((y0[tets[i].o1]-tt0)*((y1[tets[i].o3]-tt1)*(y2[tets[i].o2]-tt2) -
                   (y1[tets[i].o2]-tt1)*( y2[tets[i].o3]-tt2)) +
                   (y1[tets[i].o1]-tt1)*((y0[tets[i].o2]-tt0)*(y2[tets[i].o3]-tt2) -
                   (y0[tets[i].o3]-tt0)*( y2[tets[i].o2]-tt2)) +
                   (y2[tets[i].o1]-tt2)*((y0[tets[i].o3]-tt0)*(y1[tets[i].o2]-tt1) -
                   (y0[tets[i].o2]-tt0)*( y1[tets[i].o3]-tt1)))*sgn <1e-4)
            {
                epsilon /= 2.0;
                (*cnt) ++;
                tt0      = t0 - dh0*epsilon;
                tt1      = t1 - dh1*epsilon;
                tt2      = t2 - dh2*epsilon;
                if (++iter>32)
                {
                    (void)mexPrintf("x");
                    *hp += hp0;
                    if (!flg)
                    {
                        *hf += hf0;
                        (*n)++;
                    }
                    return;
                }
            }
        }
    }
    /* Decrease epsilon until potential is less than or equal to original potential */
    for(iter=0; iter<8; iter++)
    {
        float h1;

        tt0 = t0-dh0*epsilon;
        tt1 = t1-dh1*epsilon;
        tt2 = t2-dh2*epsilon;
        for(i=0, hp1 = 0.0; i<ntets; i++)
        {
            if (!(tets[i].mask & msk))
            {
                float Y[4][3];
                Y[0][0] = tt0*vox_f[0];            Y[0][1] = tt1*vox_f[1];            Y[0][2] = tt2*vox_f[2];
                Y[1][0] = y0[tets[i].o1]*vox_f[0]; Y[1][1] = y1[tets[i].o1]*vox_f[1]; Y[1][2] = y2[tets[i].o1]*vox_f[2];
                Y[2][0] = y0[tets[i].o2]*vox_f[0]; Y[2][1] = y1[tets[i].o2]*vox_f[1]; Y[2][2] = y2[tets[i].o2]*vox_f[2];
                Y[3][0] = y0[tets[i].o3]*vox_f[0]; Y[3][1] = y1[tets[i].o3]*vox_f[1]; Y[3][2] = y2[tets[i].o3]*vox_f[2];
                fun_ss(&(tets[i]), Y, &hp1);
            }
        }

        if (!flg)
        {
            f    = resample(fvol, tt0, tt1, tt2, dim_f);
            if (!known(f))
            {
                h0   = hp0*lambda;
                h1   = hp1*lambda;
                flg  = 1;
            }
            else
            {
                d    = f*scale-g;
                hf1  = d*d/2;
                h1   = hf1/v+hp1*lambda;
            }
        }
        else
        {
            h1   = hp1*lambda;
        }

        if (h1 <= h0)
        {
            y0[0]  = tt0;
            y1[0]  = tt1;
            y2[0]  = tt2;
            (*hp) += hp1;
            if (!flg)
            {
                *sumf += f;
                *sumg += g;
                (*hf) += hf1;
                (*n)  ++;
            }
            return;
        }
        epsilon /= 2.0;
        (*cnt) ++;
    }

    (*hp) += hp0;
    if (!flg)
    {
        *sumf += f;
        *sumg += g;
        (*hf) += hf0;
        (*n)  ++;
    }
}



static void warp3d(unsigned char g[], unsigned char f[],
    float y0[], float y1[], float y2[],
    int dim_g[3], float vox_g[3], int dim_f[3], float vox_f[3],
    float *lambda, float *epsilon, float *sigma2, float *scale, int meth)
{
    int x2, x1, x0, n, cnt,sgn;
    float hf = 0.0, hp = 0.0;
    float sumf, sumg;
    setup_consts(dim_g, vox_g);

    if ((vox_f[0]*vox_f[1]*vox_f[2]>0) == (vox_g[0]*vox_g[1]*vox_g[2]>0))
        sgn = 1;
    else
        sgn = -1;

    sumf = sumg = 0.0;
    hf = 0.0; hp = 0.0; n = 0; cnt = 0;
    for(x2=1; x2<=dim_g[2]; x2++)
    {
        if ((meth & 1) || ((x2 > 1) && (x2 < dim_g[2])))
        {
            int msk3 = 0;
            if (x2 == 1) msk3 |= 16;
            if (x2 == dim_g[2]) msk3 |= 32;

            for(x1=1; x1<=dim_g[1]; x1++)
            {
                int msk2 = msk3;
                if ((meth & 1) || ((x1 > 1) && (x1 < dim_g[1])))
                {
                    if (x1 == 1) msk2 |= 4;
                    if (x1 == dim_g[1]) msk2 |= 8;

                    if (meth & 1)
                        tweek(1, x1, x2, y0, y1, y2, dim_g, g, dim_f, vox_f, f, 1|msk2,
                            *lambda, *epsilon, *sigma2, *scale, &hf, &hp, &n, &sumf, &sumg, &cnt,sgn);
                    for(x0=2; x0<dim_g[0]; x0++)
                        tweek(x0, x1, x2, y0, y1, y2, dim_g, g, dim_f, vox_f, f, msk2,
                            *lambda, *epsilon, *sigma2, *scale, &hf, &hp, &n, &sumf, &sumg, &cnt,sgn);
                    if (meth & 1)
                        tweek(dim_g[0], x1, x2, y0, y1, y2, dim_g, g, dim_f, vox_f, f, 2|msk2,
                            *lambda, *epsilon, *sigma2, *scale, &hf, &hp, &n, &sumf, &sumg, &cnt,sgn);
                }
            }
            /* (void)mexPrintf("."); */
        }
    }
    *epsilon = *epsilon/cnt*n;
    *scale = sumg/sumf;
    *sigma2 = hf/n;
    (void)mexPrintf("%.8g\n", hf/n);

    sumf = sumg = 0.0;
    hf = 0.0; hp = 0.0; n = 0; cnt = 0;
    for(x2=dim_g[2]; x2>=1; x2--)
    {
        if ((meth & 1) || ((x2 > 1) && (x2 < dim_g[2])))
        {
            int msk3 = 0;
            if (x2 == 1) msk3 |= 16;
            if (x2 == dim_g[2]) msk3 |= 32;

            for(x1=dim_g[1]; x1>=1; x1--)
            {
                int msk2 = msk3;
                if ((meth & 1) || ((x1 > 1) && (x1 < dim_g[1])))
                {
                    if (x1 == 1) msk2 |= 4;
                    if (x1 == dim_g[1]) msk2 |= 8;

                    if (meth & 1)
                        tweek(dim_g[0], x1, x2, y0, y1, y2, dim_g, g, dim_f, vox_f, f, 2|msk2,
                            *lambda, *epsilon, *sigma2, *scale, &hf, &hp, &n, &sumf, &sumg, &cnt,sgn);
                    for(x0=dim_g[0]-1; x0>1; x0--)
                        tweek(x0, x1, x2, y0, y1, y2, dim_g, g, dim_f, vox_f, f, msk2,
                            *lambda, *epsilon, *sigma2, *scale, &hf, &hp, &n, &sumf, &sumg, &cnt,sgn);
                    if (meth & 1)
                        tweek(1, x1, x2, y0, y1, y2, dim_g, g, dim_f, vox_f, f, 1|msk2,
                            *lambda, *epsilon, *sigma2, *scale, &hf, &hp, &n, &sumf, &sumg, &cnt,sgn);
                }
            }
            /* (void)mexPrintf("."); */
        }
    }
    *epsilon = *epsilon/cnt*n;
    *scale   = sumg/sumf;
    *sigma2  = hf/n;
    (void)mexPrintf("%.8g\n", hf/n);
}

static float get_scale(unsigned char g[], unsigned char f[], float y0[], float y1[], float y2[], int dim_g[3], int dim_f[3])
{
    int x0, x1, x2, o;
    float scale = 0.0;
    float gpix, fpix;
    float sumg=0.0, sumf=0.0;
    for(x2=0; x2<dim_g[2]; x2++)
        for(x1=0; x1<dim_g[1]; x1++)
            for(x0=0; x0<dim_g[0]; x0++)
            {
                o = x0 + dim_g[0]*(x1 + x2*dim_g[1]);
                gpix = g[o];
                fpix = resample(f, y0[o], y1[o], y2[o], dim_f);
                if (known(fpix) && known(gpix))
                {
                    sumg += gpix;
                    sumf += fpix;
                }
            }

    scale = sumg/sumf;
    return(scale);
}

static float get_sumsq(unsigned char g[], unsigned char f[], float y0[], float y1[], float y2[], int dim_g[3], int dim_f[3], float scale)
{
    int x0, x1, x2, o, n=0;
    float sigma2 = 0.0, tmp;
    float gpix, fpix;

    n = 0;
    for(x2=0; x2<dim_g[2]; x2++)
        for(x1=0; x1<dim_g[1]; x1++)
            for(x0=0; x0<dim_g[0]; x0++)
            {
                o = x0 + dim_g[0]*(x1 + x2*dim_g[1]);
                gpix = g[o];
                fpix = resample(f, y0[o], y1[o], y2[o], dim_f);
                if (known(fpix) && known(gpix))
                {
                    tmp     = gpix-fpix*scale;
                    sigma2 += tmp*tmp;
                    n ++;
                }
            }
    sigma2 /= (2*n);
    return(sigma2);
}

static void estimate_warps(unsigned char g[], unsigned char f[],
    float y0[], float y1[], float y2[],
    int dim_g[3], float vox_g[3], int dim_f[3], float vox_f[3],
    int its, float lambda, float epsilon, int meth)
{
    int iter;
    float sigma2, scale;
    scale  = get_scale(g,f,y0,y1,y2,dim_g,dim_f);
    sigma2 = get_sumsq(g,f,y0,y1,y2,dim_g,dim_f, scale);
    for(iter=0; iter<its; iter++)
    {
        /* (void)mexPrintf("%d: ", iter+1); */
        warp3d(g,f,y0,y1,y2, dim_g,vox_g, dim_f,vox_f, &lambda, &epsilon, &sigma2, &scale, meth);
    }
}


static float *get_def_volume(const mxArray *ptr, int dims[3])
{
    int nd, i;
    const int *ldims;
    if (mxIsStruct(ptr) || !mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsSparse(ptr) || !mxIsSingle(ptr))
        mexErrMsgTxt("Deformations must be single precision floating point multi-dimensional arrays.");

    nd = mxGetNumberOfDimensions(ptr);
    if (nd>3)
        mexErrMsgTxt("Too many dimensions in data.");

    ldims = mxGetDimensions(ptr);
    for(i=0; i<nd; i++)
        if (dims[i] != ldims[i])
            mexErrMsgTxt("Incompatible dimensions.");
    for(i=nd; i<3; i++)
        if (dims[i] != 1)
            mexErrMsgTxt("Incompatible dimensions.");

    return((float *)mxGetPr(ptr));
}

static unsigned char *get_uint8_volume(const mxArray *ptr, int dims[3])
{
    int nd, i;
    const int *ldims;
    if (mxIsStruct(ptr) || !mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsSparse(ptr) || !mxIsUint8(ptr))
        mexErrMsgTxt("Data must be uint8 multi-dimensional arrays.");

    nd = mxGetNumberOfDimensions(ptr);
    if (nd>3)
        mexErrMsgTxt("Too many dimensions in data.");

    ldims = mxGetDimensions(ptr);
    for(i=0; i<nd; i++)
        dims[i] = ldims[i];
    for(i=nd; i<3; i++)
        dims[i] = 1;

    return((unsigned char *)mxGetPr(ptr));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned char *f=0, *g=0;
    int dim_g[3], dim_f[3], its, meth, i;
    float lambda, epsilon, vox_g[3], vox_f[3];
    float *y0, *y1, *y2;

    if (nrhs!=7 || nlhs>0) mexErrMsgTxt("Incorrect usage.");

    g = get_uint8_volume(prhs[0], dim_g);
    f = get_uint8_volume(prhs[1], dim_f);

    y0 = get_def_volume(prhs[2],dim_g);
    y1 = get_def_volume(prhs[3],dim_g);
    y2 = get_def_volume(prhs[4],dim_g);

    if (!mxIsNumeric(prhs[5]) || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
        mxGetM(prhs[5]) != 3 || mxGetN(prhs[5]) != 2)
        mexErrMsgTxt("Can't use voxel size data.");

    for(i=0; i<3; i++)
    {
        vox_g[i] = mxGetPr(prhs[5])[i];
        vox_f[i] = mxGetPr(prhs[5])[i+3];
    }

    if (!mxIsNumeric(prhs[6]) || !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
        mxGetM(prhs[6])*mxGetN(prhs[6]) != 4)
        mexErrMsgTxt("Can't use options data.");
    its     =  mxGetPr(prhs[6])[0];
    lambda  =  mxGetPr(prhs[6])[1];
    epsilon =  mxGetPr(prhs[6])[2];
    meth    =  mxGetPr(prhs[6])[3];
    estimate_warps(g,f,y0,y1,y2, dim_g, vox_g, dim_f, vox_f, its, lambda, epsilon, meth);
}
