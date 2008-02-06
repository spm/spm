/*
 * $Id: spm_def2det.c 1140 2008-02-06 19:24:05Z spm $
 * John Ashburner
 */

#include "mex.h"
#include <math.h>
#include <stdio.h>
#define REAL float

typedef struct
{
        REAL ix[4][4];
        REAL volx;
        int   o1, o2, o3;
        unsigned long mask;
} Tettype;

static Tettype oddtets[32], eventets[8];

static void invertX(REAL X[4][3], REAL IX[4][4], REAL *dt)
/* X is a matrix containing the co-ordinates of the four vertices of a tetrahedron.
   IX = inv([X ; 1 1 1 1]);  */
{
    REAL id;
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

static void mulMX(REAL A[4][3], REAL B[4][3], REAL C[4][3])
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

static void setup_consts_sub(int dim[3], int x[3][3], Tettype *tet, REAL M[4][3])
{
    REAL X[4][3], MX[4][3], dtx;
    X[0][0] = 0.0;
    X[0][1] = 0.0;
    X[0][2] = 0.0;

    X[1][0] = x[0][0];
    X[1][1] = x[0][1];
    X[1][2] = x[0][2];

    X[2][0] = x[1][0];
    X[2][1] = x[1][1];
    X[2][2] = x[1][2];

    X[3][0] = x[2][0];
    X[3][1] = x[2][1];
    X[3][2] = x[2][2];

    mulMX(MX, M, X);
    invertX(MX, tet->ix, &dtx);

    tet->volx = dtx/6.0;
    tet->o1   = x[0][0]+dim[0]*(x[0][1]+dim[1]*x[0][2]);
    tet->o2   = x[1][0]+dim[0]*(x[1][1]+dim[1]*x[1][2]);
    tet->o3   = x[2][0]+dim[0]*(x[2][1]+dim[1]*x[2][2]);
}

static void setup_consts(int dim[3], REAL M[4][3])
{
    int i;
    static int xo[][3][3] = {
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

    static int xe[][3][3] = {
        {{ 1, 0, 0}, { 0, 0, 1}, { 0, 1, 0}},
        {{-1, 0, 0}, { 0, 0, 1}, { 0,-1, 0}},
        {{-1, 0, 0}, { 0, 0,-1}, { 0, 1, 0}},
        {{ 1, 0, 0}, { 0, 0,-1}, { 0,-1, 0}},

        {{ 0, 0, 1}, {-1, 0, 0}, { 0, 1, 0}},
        {{ 0, 0, 1}, { 1, 0, 0}, { 0,-1, 0}},
        {{ 0, 0,-1}, { 1, 0, 0}, { 0, 1, 0}},
        {{ 0, 0,-1}, {-1, 0, 0}, { 0,-1, 0}}};

    for(i=0; i<32; i++) setup_consts_sub(dim, xo[i],  &oddtets[i], M);
    for(i=0; i< 8; i++) setup_consts_sub(dim, xe[i], &eventets[i], M);
}

static void get_jacmat(float *y0, float *y1, float *y2, Tettype tet, REAL J[3][3])
{
    REAL Y[4][3];
    Y[0][0] = y0[0];          Y[0][1] = y1[0];          Y[0][2] = y2[0];
    Y[1][0] = y0[tet.o1]; Y[1][1] = y1[tet.o1]; Y[1][2] = y2[tet.o1];
    Y[2][0] = y0[tet.o2]; Y[2][1] = y1[tet.o2]; Y[2][2] = y2[tet.o2];
    Y[3][0] = y0[tet.o3]; Y[3][1] = y1[tet.o3]; Y[3][2] = y2[tet.o3];

    J[0][0] = (Y[0][0]*tet.ix[0][0] + Y[1][0]*tet.ix[0][1]
             + Y[2][0]*tet.ix[0][2] + Y[3][0]*tet.ix[0][3]);
    J[0][1] = (Y[0][1]*tet.ix[0][0] + Y[1][1]*tet.ix[0][1]
             + Y[2][1]*tet.ix[0][2] + Y[3][1]*tet.ix[0][3]);
    J[0][2] = (Y[0][2]*tet.ix[0][0] + Y[1][2]*tet.ix[0][1]
             + Y[2][2]*tet.ix[0][2] + Y[3][2]*tet.ix[0][3]);

    J[1][0] = (Y[0][0]*tet.ix[1][0] + Y[1][0]*tet.ix[1][1]
             + Y[2][0]*tet.ix[1][2] + Y[3][0]*tet.ix[1][3]);
    J[1][1] = (Y[0][1]*tet.ix[1][0] + Y[1][1]*tet.ix[1][1]
             + Y[2][1]*tet.ix[1][2] + Y[3][1]*tet.ix[1][3]);
    J[1][2] = (Y[0][2]*tet.ix[1][0] + Y[1][2]*tet.ix[1][1]
             + Y[2][2]*tet.ix[1][2] + Y[3][2]*tet.ix[1][3]);

    J[2][0] = (Y[0][0]*tet.ix[2][0] + Y[1][0]*tet.ix[2][1]
             + Y[2][0]*tet.ix[2][2] + Y[3][0]*tet.ix[2][3]);
    J[2][1] = (Y[0][1]*tet.ix[2][0] + Y[1][1]*tet.ix[2][1]
             + Y[2][1]*tet.ix[2][2] + Y[3][1]*tet.ix[2][3]);
    J[2][2] = (Y[0][2]*tet.ix[2][0] + Y[1][2]*tet.ix[2][1]
             + Y[2][2]*tet.ix[2][2] + Y[3][2]*tet.ix[2][3]);

}

static float subfunk(int x0, int x1, int x2, float *Y0, float *Y1, float *Y2, int dim[3])
{
    Tettype *tets;
    int ntets, i, o;
    float dt, wt;

    if ((x0+x1+x2+1)%2){ tets =  oddtets; ntets = 32;}
    else               { tets = eventets; ntets = 8;}

    o   = x0 + dim[0]*(x1 + x2*dim[1]);

    dt = wt = 0.0;
    for(i=0; i<ntets; i++)
    {
        REAL J[3][3], d;
        get_jacmat(Y0+o, Y1+o, Y2+o, tets[i], J);

        d = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])+
            J[1][0]*(J[0][2]*J[2][1]-J[0][1]*J[2][2])+
            J[2][0]*(J[0][1]*J[1][2]-J[0][2]*J[1][1]);
        dt += tets[i].volx*d;
        wt += tets[i].volx;
    }
    return(dt/wt);
}



static void gen_dets(float y1[], float y2[], float y3[], int dim[3], float dets[])
{
    int x3, x2, x1, o;
    float NaN;
    NaN = mxGetNaN();

    for(x2=0; x2<dim[1]; x2++)
        for(x1=0; x1<dim[0]; x1++)
        {
            dets[x1 + dim[0]*(x2 + 0*dim[1])] = NaN;
        }

    for(x3=1; x3<dim[2]-1; x3++)
    {
        for(x1=0; x1<dim[0]; x1++)
            dets[x1 + dim[0]*(0 + x3*dim[1])] = NaN;

        for(x2=1; x2<dim[1]-1; x2++)
        {
            dets[0 + dim[0]*(x2 + x3*dim[1])] = NaN;
            for(x1=1; x1<dim[0]-1; x1++)
            {
                o = x1 + dim[0]*(x2 + x3*dim[1]);
                dets[o] = subfunk(x1,x2,x3,y1,y2,y3, dim);
            }
            dets[(dim[0]-1) + dim[0]*(x2 + x3*dim[1])] = NaN;

        for(x1=0; x1<dim[0]; x1++)
            dets[x1 + dim[0]*(dim[1]-1 + x3*dim[1])] = NaN;
        }
        (void)printf(".");
        (void)fflush(stdout);
    }
    for(x2=0; x2<dim[1]; x2++)
        for(x1=0; x1<dim[0]; x1++)
        {
            dets[x1 + dim[0]*(x2 + (dim[2]-1)*dim[1])] = NaN;
        }
    (void)printf("\n");
}

static float *get_volume(const mxArray *ptr, int dims[3])
{
    int nd, i;
    const int *ldims;
    if (mxIsStruct(ptr) || !mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsSparse(ptr) || !mxIsSingle(ptr))
        mexErrMsgTxt("Data must be a single precision floating point multi-dimensional array.");

    nd = mxGetNumberOfDimensions(ptr);
    if (nd>3) mexErrMsgTxt("Too many dimensions in data.");
    ldims = mxGetDimensions(ptr);
    for(i=0; i<nd; i++) dims[i] = ldims[i];
    for(i=nd; i<3; i++) dims[i] = 1;
    return((float *)mxGetPr(ptr));
}

static void get_mat(const mxArray *ptr, REAL M[4][3])
{
    int i, j;
    double *p;

        if (!mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsComplex(ptr) || !mxIsDouble(ptr) || mxGetM(ptr) != 4 || mxGetN(ptr) != 4)
                mexErrMsgTxt("Affine transform matrix must be 4x4.");
    p = (double *)mxGetPr(ptr);

    for(i=0; i<3; i++)
        for(j=0; j<4; j++)
            M[j][i] = p[i+4*j];

    if (p[3+4*0] != 0.0 || p[3+4*1] != 0.0 || p[3+4*2] != 0.0 || p[3+4*3] != 1.0)
        mexErrMsgTxt("No perspective projections allowed.");
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float  *y0=0,  *y1=0,  *y2=0, *dt=0;
    int dim[3], dimtmp[3];
    REAL M[4][3];

        if (nrhs != 4 || nlhs > 1) mexErrMsgTxt("Inappropriate usage.");

    y0 = get_volume(prhs[0], dim);
    y1 = get_volume(prhs[1], dimtmp);
    if (dim[0] != dimtmp[0] || dim[1] != dimtmp[1] || dim[2] != dimtmp[2])
        mexErrMsgTxt("Incompatible dimensions.");
    y2 = get_volume(prhs[2], dimtmp);
    if (dim[0] != dimtmp[0] || dim[1] != dimtmp[1] || dim[2] != dimtmp[2])
        mexErrMsgTxt("Incompatible dimensions.");

    get_mat(prhs[3],M);

    plhs[0] = mxCreateNumericArray(3, dim, mxSINGLE_CLASS, mxREAL);
    dt = (float *)mxGetPr(plhs[0]);

    setup_consts(dim, M);
    gen_dets(y0, y1, y2, dim, dt);
}
