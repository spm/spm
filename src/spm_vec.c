/* 
 * $Id: spm_vec.c 5720 2013-10-31 13:46:05Z guillaume $
 * Copyright 2012 Eduardo Aponte <aponteeduardo@gmail.com>
 */

#include "mex.h"
#include <math.h>
#include <string.h>

int forwardPass ( const mxArray * ica, size_t * mn, int * cf )
{
    /* Calculates the size of the output vector and checks for complex inputs
     * ica -> Input array
     * mn  -> Maximum number of elements
     * cf  -> complex flag
     */
    if ( ica == NULL ){
        return 1;
    } else if ( mxIsEmpty ( ica ) ){
        return 1;
    } else if ( mxIsCell( ica ) ) {
        int oen = 1;
        size_t i, dn;
        
        dn = mxGetNumberOfElements ( ica );
        for ( i = 0; i < dn ; i++ ){
            oen *= forwardPass( mxGetCell( ica, i), mn, cf);
        }
        return oen;
    } else if ( mxIsStruct( ica ) ){
        int i,nf,oen = 1;
        size_t dn,j;
        nf = mxGetNumberOfFields( ica );
        dn = mxGetNumberOfElements ( ica );
        for ( j = 0 ; j < nf ; j++ ){
            for ( i = 0; i < dn ; i++ ){
                oen *= forwardPass( mxGetFieldByNumber( ica, i, j), mn, cf);
            }
        }
        return oen;
    } 
    
    else if ( mxIsNumeric(ica) || mxIsLogical(ica) ){
        if ( mxIsComplex(ica) ) *cf = 1;
        *mn += mxGetNumberOfElements(ica);
        return 1;
    } else {
        return 1;
    }
}

int catNodes( const mxArray *iva , double *pr, double *pi )
{
    /* Concatenates the nodes.
     * iva -> Input vector array.
     * pr  -> Real part.
     * pi  -> Imaginary part.
     */
    size_t ts;
    int cf = (int ) mxIsComplex(iva);
    
    if ( mxIsSparse(iva) ){
        double *ipr,*ipi = NULL;
        mwIndex *ir,*jc,j;
        mwSize i;
        const mwSize *nd;    
        int t = 0;
        
        ir = mxGetIr(iva);
        jc = mxGetJc(iva);
        ipr = mxGetPr(iva);
        if ( cf ) ipi = mxGetPi (iva);
        
        nd = mxGetDimensions(iva);
        ts = (size_t ) nd[0]*nd[1];
        
        for ( i = 0; i < nd[1]; i++ ){
            for (j = jc[i]; j < jc[i+1]; j++){
                pr[i*nd[0]+ir[t]] = ipr[t];
                if ( cf ) pi[i*nd[0]+ir[t]] = ipi[t];
                t++;
            }
        }
        return nd[0]*nd[1];
    } else {
        size_t i;
        mxClassID id =  mxGetClassID( iva );
        
        ts = mxGetNumberOfElements( iva );
        
        if ( id == mxLOGICAL_CLASS ){
            mxLogical  *td;
            td = (mxLogical  *)mxGetData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double )td[i];
            }
        } else if ( id == mxINT8_CLASS ){
            signed char *td = (signed char *)mxGetData(iva);
            signed char *tdi =  NULL;
            if ( cf ) tdi = (signed char *)mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxUINT8_CLASS ){
            unsigned char * td = (unsigned char *) mxGetData(iva);
            unsigned char * tdi = NULL;
            if ( cf ) tdi = (unsigned char *) mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxUINT16_CLASS ){
            unsigned short * td = (unsigned short *)  mxGetData(iva);
            unsigned short * tdi = NULL;
            if ( cf ) tdi = (unsigned short *)  mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxINT16_CLASS ){
            short * td = (short *)  mxGetData(iva);
            short * tdi = NULL;
            if ( cf ) tdi = (short *)  mxGetImagData(iva);            
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxUINT32_CLASS ){
            unsigned int * td = (unsigned int *)  mxGetData(iva);
            unsigned int * tdi = NULL;
            if ( cf ) tdi = (unsigned int *)  mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxINT32_CLASS ){
            int * td = (int *) mxGetData(iva);
            int * tdi  = NULL;
            if ( cf ) tdi = (int *) mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxINT64_CLASS ){
            long long * td = (long long *) mxGetData(iva);
            long long * tdi = NULL;
            if ( cf ) tdi = (long long *) mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxUINT64_CLASS ){
            unsigned long long * td = (unsigned long long *) mxGetData(iva);
            unsigned long long * tdi = NULL;
            if ( cf ) tdi = (unsigned long long *) mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxSINGLE_CLASS ){
            float * td = (float *)  mxGetData(iva);
            float * tdi = NULL;
            if ( cf ) tdi = (float *)  mxGetImagData(iva);
            for (i = 0; i < ts; i++){
                pr[i]  = (double ) td[i];
                if ( cf ) pi[i]  = (double )tdi[i];
            }
        } else if ( id == mxDOUBLE_CLASS ){
            double * td = mxGetPr(iva);
            double * tdi = NULL;
            memcpy(pr,td,ts*sizeof(double));
            if ( cf ) {
                tdi = mxGetPi(iva);
                memcpy(pi,tdi,ts*sizeof(double));
            }
        }
    }
    return ts;
}

int backwardPass ( const mxArray *ica, double *pr, double *pi)
{
    if ( ica == NULL ){
        return 0;
    } else if ( mxIsEmpty(ica) ){
        return 0;
    } else if ( mxIsCell(ica) ) {
        int oen = 0;
        size_t i, dn;
        dn = mxGetNumberOfElements(ica);
        for ( i = 0; i < dn ; i++ ){
            oen += backwardPass( mxGetCell(ica,i), pr + oen , pi + oen);
        }
        return oen;
    } else if ( mxIsStruct(ica) ){
        int i,nf,oen = 0;
        size_t dn,j;
        
        nf = mxGetNumberOfFields(ica);
        dn = mxGetNumberOfElements(ica);
        
        for ( j = 0 ; j < nf ; j++ ){
            for ( i = 0; i < dn ; i++ ){
                oen += backwardPass( mxGetFieldByNumber(ica,i,j),pr + oen,pi + oen);
            }
        }
        return oen;
    } else if ( mxIsNumeric(ica) || mxIsLogical(ica) ){
        return catNodes(ica,pr,pi);
    }
    return 0;
}

void mexFunction ( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    /* Check for proper number of arguments. */
    if( nlhs > 1 ) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if( nrhs == 0 ) {
        mwSize td[2] = {0,0};
        plhs[0] = mxCreateNumericArray( 2, td, mxDOUBLE_CLASS, mxREAL);
        return;
    } else {
        size_t mn[1] = {0};
        mwSize td[2];
        double *pr,*pi;
        int di = 0;
        int ci[1] = {0};
        int ni, c;
        
        for (ni = 0; ni < nrhs; ni++){
            forwardPass(prhs[ni],mn,ci);
        }
        
        td[0] = (mwSize ) *mn;
        td[1] = 1;
        
        if ( *ci ){
            plhs[0] = mxCreateNumericArray( 2, td, mxDOUBLE_CLASS, mxCOMPLEX);
            pr      = mxGetPr(plhs[0]);
            pi      = mxGetPi(plhs[0]);
        } else {
            plhs[0] = mxCreateNumericArray( 2, td, mxDOUBLE_CLASS, mxREAL);
            pr      = mxGetPr(plhs[0]);
            pi      = NULL;
        }
        
        for ( c = 0; c < nrhs; c++){
            di += backwardPass(prhs[c], pr+di, pi+di);
        }
    }

    return;
}

