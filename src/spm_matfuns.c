/*
 * $Id: spm_matfuns.c 4452 2011-09-02 10:45:26Z guillaume $
 */

#include <math.h>

int gaussj(double *a, int n, double *b, int m)
{
    int icol = 0, ipiv[64], irow = 0, i, j, k, l, indxc[64], indxr[64], ll;
    double pivinv, big, dum;

    if (n>64) return(1);

    for (j=0; j<n; ++j)
        ipiv[j] = 0;
    for (i=0; i<n; ++i)
    {
        big = 0.0;
        for (j=0; j<n; ++j)
            if (ipiv[j] != 1)
                for (k=0; k<n; ++k)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(a[j + k * n]) >= big)
                        {
                            big  = fabs(a[j + k * n]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        return(1);
                }
        ++ipiv[icol];
        if (irow != icol)
        {
            for (l=0; l<n; ++l)
            {
                dum = a[irow + l * n];
                a[irow + l * n] = a[icol + l * n];
                a[icol + l * n] = dum;
            }
            for (l=0; l<m; ++l)
            {
                dum = b[irow + l * n];
                b[irow + l * n] = b[icol + l * n];
                b[icol + l * n] = dum;
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol + icol * n] == 0.0)
            return(1);
        pivinv = 1.0 / a[icol + icol * n];
        a[icol + icol * n] = 1.0;
        for (l=0; l<n; ++l)
            a[icol + l * n] *= pivinv;
        for (l=0; l<m; ++l)
            b[icol + l * n] *= pivinv;
        for (ll=0; ll<n; ++ll)
            if (ll != icol)
            {
                dum = a[ll + icol * n];
                a[ll + icol * n] = 0.0;
                for (l=0; l<n; ++l)
                    a[ll + l * n] -= a[icol + l * n] * dum;
                for (l=0; l<m; ++l)
                    b[ll + l * n] -= b[icol + l * n] * dum;
            }
    }
    for (l=n-1; l>=0; --l)
        if (indxr[l] != indxc[l])
            for (k=0; k<n; ++k)
            {
                dum = a[k + indxr[l] * n];
                a[k + indxr[l] * n] = a[k + indxc[l] * n];
                a[k + indxc[l] * n] = dum;
            }
    return(0);
}

/* M = MG\MF */
int AbackslashB(double *MG, double *MF, double *M)
{
    int i;
    double IMG[16];
    for(i=0; i<16; i++)
    {
        IMG[i] = MG[i];
        M[i]   = MF[i];
    }
    if (gaussj(IMG, 4, M, 4)) return(1);
    return(0);
}

void MtimesX(double *M, double *x, double *y)
{
    y[0] = M[0+4*0]*x[0] + M[0+4*1]*x[1] + M[0+4*2]*x[2] + M[0+4*3];
    y[1] = M[1+4*0]*x[0] + M[1+4*1]*x[1] + M[1+4*2]*x[2] + M[1+4*3];
    y[2] = M[2+4*0]*x[0] + M[2+4*1]*x[1] + M[2+4*2]*x[2] + M[2+4*3];
}


/*
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
{
    unsigned int n, m, i;
    double *MG, *MF, *M, tmp[16];

    if (nrhs == 0) mexErrMsgTxt("usage: M=leftdiv(MG,MF)");
    if (nrhs != 2) mexErrMsgTxt("leftdiv: 2 input arguments required");
    if (nlhs > 1) mexErrMsgTxt("leftdiv: only 1 output argument required");

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsFull(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("leftdiv: MG must be numeric, real, full and double");
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || !mxIsFull(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("leftdiv: MF must be numeric, real, full and double");

    if (mxGetM(prhs[0])!=4 || mxGetN(prhs[0])!=4 || mxGetM(prhs[1])!=4 || mxGetN(prhs[1])!=4)
        mexErrMsgTxt("leftdiv: all matrices must be 4x4");
    MG = mxGetPr(prhs[0]);
    MF = mxGetPr(prhs[1]);

    plhs[0] = mxCreateFull(4,4, REAL);
    M = mxGetPr(plhs[0]);

    AbackslashB(MG, MF, M);
}

void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
{
    unsigned int n, m, i;
    double *X, *A, *B, *IA;

    if (nrhs == 0) mexErrMsgTxt("usage: [X, IB]=leftdiv(B,C)");
    if (nrhs != 2) mexErrMsgTxt("leftdiv: 2 input arguments required");
    if (nlhs > 2) mexErrMsgTxt("leftdiv: only 2 output arguments required");

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsFull(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("leftdiv: A must be numeric, real, full and double");
    A = mxGetPr(prhs[0]);
    n = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) != n)
        mexErrMsgTxt("leftdiv: A has incompatible no of columns");
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || !mxIsFull(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("leftdiv: B must be numeric, real, full and double");
    B = mxGetPr(prhs[1]);
    if (mxGetM(prhs[1]) != n)
        mexErrMsgTxt("leftdiv: B has incompatible m dimension");
    m = mxGetN(prhs[1]);

    plhs[0] = mxCreateFull(n,m, REAL);
    X = mxGetPr(plhs[0]);

    plhs[1] = mxCreateFull(n,n, REAL);
    IA = mxGetPr(plhs[1]);

    for(i=0; i<m*n; i++) X[i]  = B[i];
    for(i=0; i<n*n; i++) IA[i] = A[i];

    (void)gaussj(IA,m,X,n);
}
*/
