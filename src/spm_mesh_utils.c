/*
 * $Id: spm_mesh_utils.c 4081 2010-10-07 14:04:44Z guillaume $
 * Guillaume Flandin
 */

#include <string.h>
#include "mex.h"

/* --- Compute geodesic distance using Dijkstra algorithm --- */
void dijkstra(double *nghbr, double *dnghbr, int nv, int nb,
              int *source, int ns,
              double maxdist,
              double *dist) {
    int i, nQ, u = 0, v, Qu = 0;
    double dmin, alt;
    int *Q = mxMalloc(nv*sizeof(int));
    
    for (i=0;i<nv;i++) {
        dist[i] = mxGetInf();
        Q[i] = i;
    }
    for (i=0;i<ns;i++)
        dist[source[i]] = 0.0;
    nQ = nv;
    
    while(nQ) {
        for (i=0,dmin=mxGetInf();i<nQ;i++) {
            v = Q[i];
            if (dist[v] < dmin) {
                dmin = dist[v];
                Qu = i;
                u = v;
            }
        }
        if ((mxIsInf(dmin)) || dmin > maxdist) break;
        Q[Qu] = Q[--nQ];
        for (i=0;i<nb;i++) {
            v = (int)nghbr[u+nv*i] - 1;
            if (v < 0) break;
            alt = dmin + dnghbr[u+nv*i];
            if (alt < dist[v])
                dist[v] = alt;
        }
    }
    mxFree(Q);
}

/* Gateway Function for Volume */
void mexFunctionVolume(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double vol, x1, x2, x3, y1, y2, y3, z1, z2, z3; 
    double *f, *v;
    mxArray *array;
    int i, nv, nf;
    
    if (nrhs == 0) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 1) mexErrMsgTxt("Too many input arguments.");
    if ((!mxIsStruct(prhs[0])) || (mxIsClass(prhs[0],"gifti")))
        mexErrMsgTxt("First argument must be a patch structure.");

    array = mxGetField(prhs[0], 0, "vertices");
    if (!mxIsDouble(array))
        mexErrMsgTxt("Vertices have to be stored as double.");
    nv    = mxGetM(array);
    v     = mxGetPr(array);
    array = mxGetField(prhs[0], 0, "faces");
    if (!mxIsDouble(array))
        mexErrMsgTxt("Faces have to be stored as double.");
    nf    = mxGetM(array);
    f     = mxGetPr(array);
    
    for (i=0,vol=0;i<nf;i++) {
        x1 = v[(int)f[i]];      y1 = v[(int)f[i]+nv];      z1 = v[(int)f[i]+2*nv];
        x2 = v[(int)f[i+nf]];   y2 = v[(int)f[i+nf]+nv];   z2 = v[(int)f[i+nf]+2*nv];
        x3 = v[(int)f[i+2*nf]]; y3 = v[(int)f[i+2*nf]+nv]; z3 = v[(int)f[i+2*nf]+2*nv];
        vol += 1.0/6.0 * (-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
    }
    vol = (vol < 0)? -vol: vol;
    plhs[0] = mxCreateDoubleScalar(vol);
}

/* Gateway Function for Neighbours */
void mexFunctionNeighbours(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize i, j, k, n, d = 0;
    mwIndex *Ir, *Jc;
    double *N = NULL, *D = NULL, *pr = NULL;
    
    if (nrhs == 0) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 1) mexErrMsgTxt("Too many input arguments.");
    if (!mxIsSparse(prhs[0])) mexErrMsgTxt("First argument must be a sparse array.");

    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), d, mxREAL);
    N = mxGetPr(plhs[0]);
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[0]), d, mxREAL);
        D = mxGetPr(plhs[1]);
    }
    
    n  = mxGetM(prhs[0]);
    Ir = mxGetIr(prhs[0]);
    Jc = mxGetJc(prhs[0]);
    pr = mxGetPr(prhs[0]);
    
    for (i=0;i<n;i++) {
        k = Jc[i+1]-Jc[i];
        if (k > d) {
            d = k; N = mxRealloc(N, mxGetM(prhs[0])*d*sizeof(double));
            if (nlhs > 1) D = mxRealloc(D, mxGetM(prhs[0])*d*sizeof(double));
        }
        for (j=0;j<k;j++) N[i+n*j] = 1 + *Ir++;
        if (nlhs > 1)
            for (j=0;j<k;j++) D[i+n*j] = *pr++;
    }
    
    mxSetPr(plhs[0],N); mxSetN(plhs[0],d);
    if (nlhs > 1) { mxSetPr(plhs[1],D); mxSetN(plhs[1],d); }
}

/* Gateway Function for Dijkstra */
void mexFunctionDijkstra(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, nv, nb, ns, *source;
    double distmax;
    
    if (nrhs < 4) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 4) mexErrMsgTxt("Too many input arguments.");
    if (!mxIsNumeric(prhs[0]))
        mexErrMsgTxt("First argument must be a neighbour array.");
    if (!mxIsNumeric(prhs[1]))
        mexErrMsgTxt("Second argument must be a distance array.");
    
    nv = mxGetM(prhs[0]);
    nb = mxGetN(prhs[0]);
    
    ns = mxGetNumberOfElements(prhs[2]);
    source = mxMalloc(ns*sizeof(double));
    for (i=0;i<ns;i++) {
        source[i] = (int)mxGetPr(prhs[2])[i];
        if ((source[i]<0) || (source[i]>=nv)) mexErrMsgTxt("Invalid vertex index.");
    }
    
    distmax = mxGetScalar(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(nv, 1, mxREAL);
    
    dijkstra(mxGetPr(prhs[0]), mxGetPr(prhs[1]), nv, nb, source, ns, distmax, mxGetPr(plhs[0]));
    
    mxFree(source);
}

/* Main Gateway Function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *action = NULL;
    mwSize strlgh;
    
    if (nrhs < 2) mexErrMsgTxt("Not enough input arguments.");
    if (!mxIsChar(prhs[0])) mexErrMsgTxt("First argument must be an action string.");
    
    strlgh = (mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(char)) + 1;
    action = mxCalloc(strlgh, sizeof(char));
    mxGetString(prhs[0], action, strlgh);
    
    if (!strcmp(action,"dijkstra")) {
        mexFunctionDijkstra(nlhs,plhs,nrhs-1,&prhs[1]);
    }
    else if (!strcmp(action,"neighbours")) {
        mexFunctionNeighbours(nlhs,plhs,nrhs-1,&prhs[1]);
    }
    else if (!strcmp(action,"volume")) {
        mexFunctionVolume(nlhs,plhs,nrhs-1,&prhs[1]);    
    }
    else {
        mexErrMsgTxt("Unknown action.");
    }
    
    mxFree(action);
}
