/*
 * $Id: spm_mesh_reduce.cpp 7421 2018-09-20 10:58:01Z guillaume $
 * Guillaume Flandin
 */

#include "mex.h"
#include "Simplify.h"

/* Gateway Function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *array = NULL;
    void *f = NULL, *v = NULL;
    double *out = NULL;
    unsigned int i;
    double agressiveness = 7.0;
    int target = 0;
    mwSize nv, nf;
    bool isVdouble = true, isFdouble = true;
    const char *fnames[] = {"vertices", "faces"};
    
    if (nrhs == 0) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 3) mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    
    if ((!mxIsStruct(prhs[0])) || (mxIsClass(prhs[0],"gifti")))
        mexErrMsgTxt("First argument must be a patch structure.");
    
    if (nrhs > 1) target = (int)mxGetScalar(prhs[1]);
    
    if (nrhs > 2) agressiveness = mxGetScalar(prhs[2]);
    
    array = mxGetField(prhs[0], 0, "vertices");
    if (array == NULL)
        mexErrMsgTxt("Field 'vertices' missing.");
    else if (!mxIsDouble(array) && !mxIsClass(array,"single"))
        mexErrMsgTxt("Vertices have to be stored as floating point numbers.");
    nv    = mxGetM(array);
    v     = mxGetData(array);
    
    array = mxGetField(prhs[0], 0, "faces");
    if (array == NULL)
        mexErrMsgTxt("Field 'faces' missing.");
    else if (!mxIsDouble(array) && !mxIsClass(array,"int32"))
        mexErrMsgTxt("Faces have to be stored as double or int32.");
    nf    = mxGetM(array);
    f     = mxGetData(array);
    
    Simplify::vertices.clear();
    Simplify::Vertex vert;
    isVdouble = mxIsDouble(mxGetField(prhs[0], 0, "vertices"));
    for (i=0;i<nv;i++) {
        if (isVdouble) {
            vert.p.x = ((double *)v)[i];
            vert.p.y = ((double *)v)[i+nv];
            vert.p.z = ((double *)v)[i+2*nv];
        } else {
            vert.p.x = ((float *)v)[i];
            vert.p.y = ((float *)v)[i+nv];
            vert.p.z = ((float *)v)[i+2*nv];
        }
        Simplify::vertices.push_back(vert);
    }
    
    Simplify::triangles.clear();
    Simplify::Triangle t;
    isFdouble = mxIsDouble(mxGetField(prhs[0], 0, "faces"));
    t.attr = 0;
    t.material = -1;
    for (i=0;i<nf;i++) {
        if (isFdouble) {
            t.v[0] = ((double *)f)[i]-1;
            t.v[1] = ((double *)f)[i+nf]-1;
            t.v[2] = ((double *)f)[i+2*nf]-1;
        } else {
            t.v[0] = ((int *)f)[i]-1;
            t.v[1] = ((int *)f)[i+nf]-1;
            t.v[2] = ((int *)f)[i+2*nf]-1;
        }
        Simplify::triangles.push_back(t);
    }

    if (target == 0) target = (int) (Simplify::triangles.size() / 2);
    
    Simplify::simplify_mesh(target, agressiveness, false);
    
    plhs[0] = mxCreateStructMatrix(1, 1, 2, fnames);
    array = mxCreateNumericMatrix(Simplify::vertices.size(), 3, mxDOUBLE_CLASS, mxREAL);
    out = mxGetPr(array);
    for (i=0;i<Simplify::vertices.size();i++) {
        out[i] = Simplify::vertices[i].p.x;
        out[i+Simplify::vertices.size()] = Simplify::vertices[i].p.y;
        out[i+2*Simplify::vertices.size()] = Simplify::vertices[i].p.z;
    }
    mxSetFieldByNumber(plhs[0], 0, 0, array);
    
    array = mxCreateNumericMatrix(Simplify::triangles.size(), 3, mxDOUBLE_CLASS, mxREAL);
    out = mxGetPr(array);
    for (i=0;i<Simplify::triangles.size();i++) {
        if (!Simplify::triangles[i].deleted) {
            out[i] = Simplify::triangles[i].v[0]+1;
            out[i+Simplify::triangles.size()] = Simplify::triangles[i].v[1]+1;
            out[i+2*Simplify::triangles.size()] = Simplify::triangles[i].v[2]+1;
        }
    }
    mxSetFieldByNumber(plhs[0], 0, 1, array);
}
