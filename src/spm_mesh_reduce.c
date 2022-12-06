/*
 * Guillaume Flandin
 * Copyright (C) 2022 Wellcome Centre for Human Neuroimaging
 */

#include "mex.h"
#include "external/nii2mesh/quadric.h"

/* Gateway Function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *array = NULL;
    void *f = NULL, *v = NULL;
    double *out = NULL;
    int i, nv, nf, target = 0;
    double aggressiveness = 7.0;
    bool isVdouble = true, isFdouble = true;
    const char *fnames[] = {"vertices", "faces"};
    vec3d *vertices = NULL;
	vec3i *faces = NULL;
    
    if (nrhs == 0) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 3) mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    
    if ((!mxIsStruct(prhs[0])) || (mxIsClass(prhs[0],"gifti")))
        mexErrMsgTxt("First argument must be a patch structure.");
    
    if (nrhs > 1) target = (int)mxGetScalar(prhs[1]);
    
    if (nrhs > 2) aggressiveness = mxGetScalar(prhs[2]);
    
    array = mxGetField(prhs[0], 0, "vertices");
    if (array == NULL)
        mexErrMsgTxt("Field 'vertices' missing.");
    else if (!mxIsDouble(array) && !mxIsClass(array,"single"))
        mexErrMsgTxt("Vertices have to be stored as floating point numbers.");
    nv    = (int)mxGetM(array);
    v     = mxGetData(array);
    
    array = mxGetField(prhs[0], 0, "faces");
    if (array == NULL)
        mexErrMsgTxt("Field 'faces' missing.");
    else if (!mxIsDouble(array) && !mxIsClass(array,"int32"))
        mexErrMsgTxt("Faces have to be stored as double or int32.");
    nf    = (int)mxGetM(array);
    f     = mxGetData(array);
    
    isVdouble = mxIsDouble(mxGetField(prhs[0], 0, "vertices"));
    vertices = (vec3d*)malloc(nv * sizeof(vec3d));
    for (i=0;i<nv;i++) {
        if (isVdouble) {
            vertices[i].x = ((double *)v)[i];
            vertices[i].y = ((double *)v)[i+nv];
            vertices[i].z = ((double *)v)[i+2*nv];
        } else {
            vertices[i].x = ((float *)v)[i];
            vertices[i].y = ((float *)v)[i+nv];
            vertices[i].z = ((float *)v)[i+2*nv];
        }
    }

    isFdouble = mxIsDouble(mxGetField(prhs[0], 0, "faces"));
    faces = (vec3i*)malloc(nf * sizeof(vec3i));
    for (i=0;i<nf;i++) {
        if (isFdouble) {
            faces[i].x = ((double *)f)[i]-1;
            faces[i].y = ((double *)f)[i+nf]-1;
            faces[i].z = ((double *)f)[i+2*nf]-1;
        } else {
            faces[i].x = ((int *)f)[i]-1;
            faces[i].y = ((int *)f)[i+nf]-1;
            faces[i].z = ((int *)f)[i+2*nf]-1;
        }
    }
    
    if (target == 0) target = (int) (nf / 2);

    quadric_simplify_mesh(&vertices, &faces, &nv, &nf, target, aggressiveness, false, aggressiveness < 7.0);

    plhs[0] = mxCreateStructMatrix(1, 1, 2, fnames);
    
    array = mxCreateNumericMatrix(nv, 3, mxDOUBLE_CLASS, mxREAL);
    out = mxGetPr(array);
    for (i=0;i<nv;i++) {
        out[i]      = vertices[i].x;
        out[i+nv]   = vertices[i].y;
        out[i+2*nv] = vertices[i].z;
    }
    mxSetFieldByNumber(plhs[0], 0, 0, array);
    
    array = mxCreateNumericMatrix(nf, 3, mxDOUBLE_CLASS, mxREAL);
    out = mxGetPr(array);
    for (i=0;i<nf;i++) {
        out[i]      = faces[i].x + 1;
        out[i+nf]   = faces[i].y + 1;
        out[i+2*nf] = faces[i].z + 1;
    }
    mxSetFieldByNumber(plhs[0], 0, 1, array);
    
    free(vertices);
    free(faces);
}
