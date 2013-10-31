/* 
 * $Id: spm_unvec.c 5720 2013-10-31 13:46:05Z guillaume $
 * Copyright 2012 Eduardo Aponte <aponteeduardo@gmail.com>
 */

#include "mex.h"
#include <math.h>
#include <string.h>

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

/* Full to sparse */
int full2sparse( const mxArray *iva, const mxArray *ia, mxArray **oa, int tt )
{
    /* Convert an array into a sparse matrix */
    mwIndex i, j, k = 0;
    double *opr, *opi, *pr, *pi, percent_sparse = 0.5;
    mwSize nzmax, *ojc, *oir;
    const mwSize *ds;
    int cf = (int ) mxIsComplex( iva );
    size_t ne = mxGetNumberOfElements( ia );

    /* Get the size and pointers to input data */
    
    ds = mxGetDimensions( ia );
    
    pr = mxGetPr(iva) + tt;
    if ( cf ) pi = mxGetPi(iva) + tt ;
    else pi = NULL;
    
    if ( ne == 0 ){
        if ( cf ) *oa = mxCreateSparse(ds[0], ds[1], 0, mxCOMPLEX);
        else *oa = mxCreateSparse(ds[0], ds[1], 0, mxREAL);
        return 0;
    }

    
    nzmax = (mwSize ) ceil((double) ne * percent_sparse);
    opr  = mxCalloc(nzmax,sizeof(double));
    ojc  = mxMalloc((ds[1]+1)*sizeof(mwIndex));
    oir  = mxMalloc(nzmax*sizeof(mwIndex));
    
    if ( cf ) opi = mxCalloc(nzmax,sizeof(double));
    else opi = NULL;
    
    /* Copy nonzeros */
    
    for ( j=0 ; ( j < ds[1] ); j++ ) {
        ojc[j] = k;
        for ( i = 0 ; i < ds[0]; i++) {
            int t = 0;
            if ( IsNonZero( pr[i+(j*ds[0])] ) ){
                t = 1;
            } else if ( cf ) {
                if ( IsNonZero( pi[i+(j*ds[0])] ) ) t = 1;
            }
           
            if ( t ) {
                /* Check to see if non-zero element will fit in
                 * allocated output array.  If not, increase percent_sparse
                 * by 10%, recalculate nzmax, and augment the sparse array */
                
                if ( k >= nzmax){
                    mwSize oldnzmax = nzmax;
                    percent_sparse += 0.5;
                    nzmax = (mwSize)ceil( (double)ds[0]*(double)ds[1]*percent_sparse);
                    
                    /* Make sure nzmax increases at least by 1 */
                    if (oldnzmax == nzmax){
                        nzmax++;
                    }
                    opr =  mxRealloc(opr, nzmax*sizeof(double));
                    oir =  mxRealloc(oir, nzmax*sizeof(mwSize));
                    if ( cf ) opi = mxRealloc(opi, nzmax*sizeof(double));
                }
                opr[k] = pr[i+(j*ds[0])];
                oir[k] = (mwSize )i;
                if ( cf ) opi[k] = pi[i+(j*ds[0])];
                k++;
            }
        }
    }

    
    ojc[ds[1]] = (mwSize ) k;
    
    if ( cf ) *oa = mxCreateSparse(ds[0],ds[1],(mwSize ) k,mxCOMPLEX);
    else *oa = mxCreateSparse(ds[0],ds[1],(mwSize ) k,mxREAL);
    
    if ( cf ){
        mxFree(mxGetPi(*oa));
        mxSetPi(*oa,opi);
    }    
    mxFree(mxGetPr(*oa));
    mxSetPr(*oa,opr);
    mxFree(mxGetJc(*oa));
    mxSetJc(*oa,ojc);
    mxFree(mxGetIr(*oa));
    mxSetIr(*oa,oir);
    
    return tt + ne;
}

int enterNode( mxArray *iva, const mxArray * ica, mxArray ** oca, int t)
{    
    if (ica == NULL){
        *oca = NULL;
        return t;
    } else if ( mxIsCell(ica) ){
        size_t i;
        const mwSize * d = mxGetDimensions(ica);
        mwSize dn = mxGetNumberOfDimensions(ica);
        size_t ne = mxGetNumberOfElements(ica);

        *oca = mxCreateCellArray(dn,d);
        
        for ( i = 0; i < ne ; i++ ){
            mxArray * toca;
            t = enterNode( iva, mxGetCell(ica,i), &toca , t );
            mxSetCell(*oca,i,toca);
        }
        return t;
    } else if ( mxIsStruct(ica) ){
        const mwSize *d = mxGetDimensions(ica);
        const size_t ne = mxGetNumberOfElements (ica);
        const char **names;
        mwSize dn = mxGetNumberOfDimensions(ica);
        mwIndex i;
        int fn = mxGetNumberOfFields(ica);
        int j;
        
        names = mxMalloc ( fn * sizeof(char*));

        for ( i = 0; i < (mwIndex ) fn ; i ++){
            names[i] = mxGetFieldNameByNumber(ica,(int)i);
        }
        *oca = mxCreateStructArray(dn,d,fn,names);
        for ( j = 0 ; j < fn ; j++ ){
            for ( i = 0; i < ne ; i++ ){
                mxArray *toca;
                t = enterNode(iva,mxGetFieldByNumber(ica,i,j),&toca,t);
                mxSetFieldByNumber(*oca,i,j,toca);
            }
        }
        return t;
    } else if ( mxIsSparse(ica) ){
        size_t ne = mxGetNumberOfElements(ica);
        if ( (ne + (size_t ) t) > mxGetNumberOfElements(iva) ){
            mexErrMsgTxt("Number of elements in vector is too high.");
        }
        full2sparse(iva,ica,oca,t);
        return t + ne;
    } else if ( mxIsDouble(ica) ){
        mwSize dn = mxGetNumberOfDimensions(ica);
        const mwSize *d = mxGetDimensions(ica);
        mxComplexity cf;
        double *pr, *pi, *tpi;
        int ne = mxGetNumberOfElements(ica);
        
        if ( ne + t > mxGetNumberOfElements(iva) ){
            mexErrMsgTxt("Number of elements in vector is too high.");
        }
        
        if ( mxIsComplex(iva) ) cf = mxCOMPLEX;
        else cf = mxREAL;
        
        *oca = mxCreateNumericArray(dn,d,mxDOUBLE_CLASS,cf );
        
        tpi = (double *) mxGetData(iva);
        
        mxFree(mxGetPr(*oca));
        pr = mxCalloc(ne,sizeof(double));
        memcpy(pr,tpi + t,ne*sizeof(double));
        mxSetPr(*oca,pr);
        
        if ( mxIsComplex(ica) ){
            mxFree(mxGetPi(*oca));
            pi = mxMalloc(ne*sizeof(double));
            memcpy(pi,tpi + t,ne*sizeof(double));
            mxSetPi(*oca,pi);
        }
        return t + ne;
    } else {
        mwSize d[2] = {0,0};
        *oca = mxCreateNumericArray(2,d,mxDOUBLE_CLASS,mxREAL);
        return t;        
    }   
}

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    mxArray *vX;
    int t;

    /* Check for proper number of arguments */
    if ( nrhs < 1 ){
        mexErrMsgTxt("vX is not defined");
    }
    if ( nrhs - 1 < nlhs ){
        mexErrMsgTxt("Too many output arguments.");        
    }
    
    t = mexCallMATLAB(1, &vX, 1, (mxArray **) prhs, "spm_vec");
    if ( t ){
        mexErrMsgTxt("Error while vectorizing the input argument");
    }
        
    if ( nlhs == 1 && nrhs == 2 ) {
        enterNode(vX,prhs[1],plhs,0);
    } else if (nlhs == 1 && nrhs > 2) {
        int t = 0;
        mwIndex i;
        mwSize td[2] = {1,1};
        td[1] = (mwSize ) nrhs-1;
        
        plhs[0] = mxCreateCellArray(2,td);
        
        for ( i = 0; i < (mwIndex ) nrhs - 1; i++ ){
            mxArray *tcoa;
            t = enterNode(vX,prhs[i+1],&tcoa,t);
            mxSetCell(*plhs,i,tcoa);
        }
    } else {
        int t = 0;
        mwIndex i;
        mwSize td[2] = {1};
        
        td[1] = (mwSize ) nrhs-1;
        for ( i = 0; i < (mwIndex ) nlhs; i++ ){
            t = enterNode(vX,prhs[i+1],plhs + i,t);
        }
    }
    
    mxDestroyArray(vX);
        
    return;
}

