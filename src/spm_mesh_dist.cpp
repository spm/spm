/*
 * Guillaume Flandin
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include <array>
#include <vector>

#include "mex.h"
#include "external/TriangleMeshDistance/TriangleMeshDistance.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *array = NULL;
    double *f = NULL, *v = NULL, *xyz = NULL;
    mwSize nv, nf, n;
    mwIndex i;
    bool signed_dist = true;
    
    /* Check for proper number of arguments. */
    if (nrhs < 2) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 3) mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    
    if ((!mxIsStruct(prhs[0])) || (mxIsClass(prhs[0],"gifti")))
        mexErrMsgTxt("First argument must be a patch structure.");
    
    /* Get vertex list */
    array = mxGetField(prhs[0], 0, "vertices");
    if (!array)
        mexErrMsgTxt("Field 'vertices' missing from patch structure.");
    if (!mxIsDouble(array) || mxGetN(array) != 3)
        mexErrMsgTxt("Vertices have to be stored as an nv x 3 double array.");
    nv    = mxGetM(array);
    v     = mxGetPr(array);
    
    /* Get triangle list */
    array = mxGetField(prhs[0], 0, "faces");
    if (!array)
        mexErrMsgTxt("Field 'faces' missing from patch structure.");
    if (!mxIsDouble(array) || mxGetN(array) != 3)
        mexErrMsgTxt("Faces have to be stored as an nf x 3 double array.");
    nf    = mxGetM(array);
    f     = mxGetPr(array);

    /* Load vertex list */
    std::vector<std::array<double, 3>> vertices;
    for (i=0;i<nv;i++)
        vertices.push_back({ v[i], v[i+nv], v[i+2*nv] });
    
    /* Load face list */
    std::vector<std::array<int, 3>> faces;
    for (i=0;i<nf;i++)
        faces.push_back({ (int)f[i]-1, (int)f[i+nf]-1, (int)f[i+2*nf]-1 });
    
    /* Construct a new TriangleMeshDistance object */
    tmd::TriangleMeshDistance mesh_distance(vertices, faces);
    
    /* Get coordinate list */
    if (!mxIsDouble(prhs[1]) || mxGetN(prhs[1]) != 3)
        mexErrMsgTxt("XYZ coordinates have to be stored as an n x 3 double array.");
    xyz = mxGetPr(prhs[1]);
    n = mxGetM(prhs[1]);

    /* Get signed/unsigned */
    if (nrhs > 2)
    {
        signed_dist = mxIsLogicalScalarTrue(prhs[2]);
    }

    /* Computes the signed distance from coordinates to the triangle mesh */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    
    tmd::Result result;
    for (i=0;i<n;i++)
    {
        if (signed_dist)
            result = mesh_distance.signed_distance({ xyz[i], xyz[i+n], xyz[i+2*n] });
        else
            result = mesh_distance.unsigned_distance({ xyz[i], xyz[i+n], xyz[i+2*n] });
        mxGetPr(plhs[0])[i] = result.distance;
    }
        
}
