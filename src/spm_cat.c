/*
 * $Id: spm_cat.c 7532 2019-02-14 12:03:24Z guillaume $
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

struct s_m {
    mwSize *jc;
    mwSize *ir;
    double *pr;
    double *pi;
    mwSize nc;
    mwSize nr;
    mwSize nz;
    int EMPTY;
    int COMP;
};


int concatenate( const mxArray *ia, struct s_m *osm );

void cast2double( const mxArray *iva, double *pr, double *pi){
    size_t i;
    size_t ts    = mxGetNumberOfElements(iva);
    mxClassID id = mxGetClassID(iva);
    int fi       = (int ) mxIsComplex(iva);
    
    if ( id == mxLOGICAL_CLASS ){
        mxLogical  *td;
        td = (mxLogical  *)mxGetData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
        }
        return;
    } else if ( id == mxINT8_CLASS ){
        signed char *tdi = NULL;
        signed char *td = (signed char * ) mxGetData(iva);
        if ( fi ) tdi = (signed char * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    } else if ( id == mxUINT8_CLASS ){
        unsigned char *td = (unsigned char * ) mxGetData(iva);
        unsigned char *tdi = NULL;
        if ( fi ) tdi = (unsigned char * ) mxGetImagData(iva);
        else tdi = NULL;
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    } else if ( id == mxUINT16_CLASS ){
        unsigned short *td = (unsigned short * ) mxGetData(iva);
        unsigned short *tdi = NULL;
        if ( fi ) tdi = (unsigned short * ) mxGetImagData(iva);

        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    } else if ( id == mxINT16_CLASS ){
        short *td = (short * ) mxGetData(iva);
        short *tdi = NULL;
        if ( fi ) tdi = (short * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double ) tdi[i];
        }
        return;        
    } else if ( id == mxUINT32_CLASS ){
        unsigned int *td = (unsigned int * ) mxGetData(iva);
        unsigned int *tdi = NULL;
        if ( fi ) tdi = (unsigned int * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double ) tdi[i];
        }
        return;
    } else if ( id == mxINT32_CLASS ){
        int *td = (int * ) mxGetData(iva);
        int *tdi = NULL;
        if ( fi ) tdi = (int * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    } else if ( id == mxINT64_CLASS ){
        long long *td = (long long * ) mxGetData(iva);
        long long *tdi = NULL;
        if ( fi ) tdi = (long long * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    } else if ( id == mxUINT64_CLASS ){
        unsigned long long * td = (unsigned long long * ) mxGetData(iva);
        unsigned long long * tdi = NULL;
        if ( fi ) tdi = (unsigned long long * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    } else if ( id == mxSINGLE_CLASS ){
        float * td = (float * ) mxGetData(iva);
        float * tdi = NULL;
        if ( fi ) tdi = (float * ) mxGetImagData(iva);
        for (i = 0; i < ts; i++){
            pr[i]  = (double ) td[i];
            if ( fi ) pi[i]  = (double )tdi[i];
        }
        return;
    }
    mexErrMsgTxt("Input is not supported. Error ocurrer while casting input into double.");
}

void free_s_m( struct s_m ism )
{
    mxFree(ism.jc);
    mxFree(ism.ir);
    mxFree(ism.pr);
    if ( ism.COMP ) mxFree(ism.pi);
}

void e_s_m( struct s_m *ism )
{
    /* Starts an empty s_m structure */
    
    ism->nc   = 0;
    ism->nr   = 0;
    ism->nz   = 0;
    /* These fields are alway allocated to prevent free empty structures */
    ism->jc = mxCalloc(1,sizeof(mwSize));
    ism->ir = mxCalloc(1,sizeof(mwSize));
    ism->pr = mxCalloc(1,sizeof(double));
    ism->EMPTY = 1;
    ism->COMP  = 0;
}

void tmc_s_m( const int cf, const mwSize nr, 
        mwSize *ijc, mwSize *iir, double *ipr, double *ipi,
        mwSize *ojc, mwSize *oir, double *opr, double *opi)
{
    /* Transfers the memory contents from one sparse matrix to other.
     * cf       Flag for complex matrices.
     * nr       Number of Rows
     * ijc      Input Jc
     * iir      Input Ir
     * ipr      Input real part
     * ipi      Input imaginary part
     * ojc      Output Jc
     * oir      Output Ir
     * opr      Output real part
     * opi      Output imaginary part */
    
    memcpy(ojc,ijc,(nr+1)*sizeof(mwSize));
    memcpy(oir,iir,ijc[nr]*sizeof(mwSize));
    memcpy(opr,ipr,ijc[nr]*sizeof(double));
    if ( cf ) memcpy(opi,ipi,ijc[nr]*sizeof(double));
}


void print_s_m( struct s_m ism )
{
    /* For debugging */
    mwSize i;
    mexPrintf("\n");
    mexPrintf("Is complex       : %4d\n",ism.COMP);
    mexPrintf("Number of columns: %4d\n",ism.nc);
    mexPrintf("Number of rows:    %4d\n",ism.nr);
    mexPrintf("Number of naz nz: %4d\n",ism.nz);
    mexPrintf("Number of nnz:     %4d\n",ism.jc[ism.nr]);
    mexPrintf("--------------------------------\n");
    for (  i = 0; i < ism.nr+1; i++ ){
        mexPrintf("%4d ",ism.jc[i]);
    }
    mexPrintf("\n");
    for (  i = 0; i < ism.nr+1; i++ ){
        mexPrintf(" --- ");
    }
    mexPrintf("\n");
    if ( ism.jc[ism.nr] > 0 )
    {
        for (  i = 0; i < ism.nr; i++ ){
            mexPrintf("%4d ",ism.ir[ism.jc[i]]);
        }
        mexPrintf("\n");
        for (  i = 0; i < ism.nr+1; i++ ){
            mexPrintf(" --- ");
        }
    }
}

void sparse2s_m( const mxArray *mxs , struct s_m *osm )
{
    mwSize *jc, *ir;
    const mwSize *ds;
    double *pr, *pi = NULL;
    
    ds = mxGetDimensions(mxs);
        
    pr = mxGetPr(mxs);
    jc = mxGetJc(mxs);
    ir = mxGetIr(mxs);
    
    osm->nc = ds[0];
    osm->nr = ds[1];
    
    osm->nz = jc[ds[1]];
    
    osm->pr = mxCalloc(osm->nz,sizeof(double));
    osm->ir = mxCalloc(osm->nz,sizeof(mwSize));
    osm->jc = mxCalloc(osm->nr+1,sizeof(mwSize));
    osm->COMP = 0;
    
    if ( mxIsComplex( mxs ) ) {
        pi = mxGetPi( mxs );
        osm->COMP = 1;
        osm->pi  = mxCalloc(osm->nz,sizeof(double));
    }
    
    tmc_s_m((int ) mxIsComplex(mxs),ds[1],
            jc,ir,pr,pi,osm->jc,osm->ir,osm->pr,osm->pi);
    
}
/* Full to sparse */
void full2sparse( const mxArray *ia, struct s_m *osm )
{
    /* Convert an array into a sparse matrix */
    
    mwSize nzmax;
    const mwSize *ds;
    mwIndex i, j, k;
    double *pr, *pi = NULL, percent_sparse;
    int cf = (int ) mxIsComplex( ia );
    
    ds = mxGetDimensions( ia );
    
    if ( mxIsDouble ( ia ) ){
        pr = mxGetPr(ia);
        if ( cf ) pi = mxGetPi(ia);
    } else {
        size_t ts = mxGetNumberOfElements( ia );
        pr = mxCalloc(ts,sizeof(double));
        if ( cf ) pi = mxCalloc(ts,sizeof(double));
        cast2double(ia,pr,pi);
    }
        
    percent_sparse = 0.5;
    
    nzmax = (mwSize ) ceil((double ) ds[1] * (double ) ds[0] * percent_sparse);
    
    osm->pr  = mxCalloc(nzmax,sizeof(double));
    osm->jc  = mxCalloc(ds[1]+1,sizeof(mwIndex));
    osm->ir  = mxCalloc(nzmax,sizeof(mwIndex));
    if ( cf ) osm->pi = mxCalloc(nzmax,sizeof(double));
    
    /* Copy nonzeros */
    
    k      = 0;
    
    for ( j=0 ; ( j < ds[1] ); j++ ) {
        ( osm -> jc )[j] = k;
        for ( i = 0 ; i < ds[0]; i++) {
            bool t = 0;
            if ( IsNonZero( pr[i+(j*ds[0])] ) ){
                t = 1;
            } else if ( cf ) {
                if ( IsNonZero( pi[i+(j*ds[0])] ) ){
                    t = 1;
                }
            }
            
            if ( t ) {
                /* Check to see if non-zero element will fit in
                 * allocated output array.  If not, increase percent_sparse
                 * by 10%, recalculate nzmax, and augment the sparse array */
                if ( k >= nzmax){
                    mwSize oldnzmax = nzmax;
                    percent_sparse += 0.5;
                    nzmax = (mwSize)ceil( (double)ds[0]*(double)ds[1]
                            *percent_sparse);
                    
                    /*make sure nzmax increases atleast by 1*/
                    if (oldnzmax == nzmax){
                        nzmax++;
                    }
                    osm->pr =  mxRealloc(osm->pr,nzmax*sizeof(double));
                    osm->ir =  mxRealloc(osm->ir,nzmax*sizeof(mwSize));
                    if ( cf ) osm->pi = mxRealloc(osm->pi,nzmax*sizeof(double));
                }
                osm->pr[k] = pr[i+(j*ds[0])];
                osm->ir[k] = (mwSize ) i;
                if ( cf ) osm->pi[k] = pi[i+(j*ds[0])];
                k++;
            }
        }
    }
    
    osm->jc[ds[1]] = (mwSize ) k;
    osm->nc = ds[0];
    osm->nr = ds[1];
    osm->nz = (mwSize ) k;
    osm->COMP = cf;
    
    /* If necessary free the memory used for casting types */
    if ( !mxIsDouble ( ia ) ){
        mxFree( pr );
        if ( cf ) mxFree( pi );
    }
    
    return;
}

mwSize cumsum(mwSize *ia,const mwIndex mt)
{
    mwSize oi = 0;
    mwIndex i;
    for (i = 0; i < mt; i++){
        oi += ia[i];
    }
    return oi;
}

mwSize forwardPass( const mxArray *ia, 
        struct s_m *tcia, mwSize *dc, mwSize *dr )
{    
    mwIndex i, j;
    mwSize dn, nz = 0;
    const mwSize *ds;
    
    ds = mxGetDimensions(ia);
    dn = mxGetNumberOfDimensions(ia);

    if ( dn > 2 ) {
        mexErrMsgTxt("Only two dimensional arrays are supported.");
    }
        
    for (i = 0; i < ds[0]; i++){
        for (j = 0; j < ds[1] ; j++){
            mwIndex ti   = i + j*ds[0];
            mxArray *tia = mxGetCell(ia,ti);
            
            if (tia == NULL){
                e_s_m( tcia + ti);
            } else if ( mxIsCell(tia) ){
                tcia[ti].EMPTY = 0;
                tcia[ti].COMP  = 0;
                concatenate( tia , tcia+ti );
            } else if  ( mxIsSparse(tia) ){
                if ( mxGetNumberOfElements(tia) ){
                    tcia[ti].EMPTY = 0;
                } else {
                    tcia[ti].EMPTY = 1;
                }
                tcia[ti].COMP  = 0;
                sparse2s_m( tia , tcia+ti );
                
            } else if ( mxIsNumeric(tia) || mxIsLogical(tia)) {
                if ( mxGetNumberOfElements(tia) ){
                    tcia[ti].EMPTY = (mwSize ) 0;
                } else {
                    tcia[ti].EMPTY = (mwSize ) 1;
                }
                tcia[ti].COMP  = 0;
                full2sparse(tia,tcia+ti);
            } else {
                mexErrMsgTxt("Only numeric and logical types are supported.");                                
            }
            
            /* Checks that the internal arrays are compatible */
            
            if ( tcia[ti].nc == 0 ){
                ;
            } else if ( dc[i+1] == 0 ){
                dc[i+1] = tcia[ti].nc;
            } else if ( tcia[ti].nc != dc[i+1]){
                mexErrMsgTxt("Matrices don't match.");
            }
            
            if ( tcia[ti].nr == 0 ){
                ;
            } else if ( dr[j+1] == 0 ){
                dr[j+1] = tcia[ti].nr;
            } else if ( tcia[ti].nr != dr[j+1]){
                mexErrMsgTxt("Matrix don't match.");
            }
            nz += tcia[ti].nz;
        }
    }
    
    return nz;
}

void backwardPass( const struct s_m *tcia, struct s_m *osm, 
        mwSize *dc, mwSize *dr, const mwSize *ds)
{
    mwSize i, j, nr;
    mwIndex toi = 0;
    
    for ( j = 0 ; j < ds[1] ; j++ ){
        for ( nr = 0; nr < dr[j+1]; nr++){
            int ts = 1;
            for ( i = 0; i < ds[0]; i++ ){
                mwIndex nl, ti = i + j*ds[0];
                /* Case of empty array */
                if ( tcia[ti].EMPTY ){
                    continue ;
                }
                
                for ( nl = tcia[ti].jc[nr] ; nl < tcia[ti].jc[nr+1] ; nl++ ){
                    osm->pr[toi] = tcia[ti].pr[nl];
                    if ( tcia[ti].COMP ) osm->pi[toi] = tcia[ti].pi[nl];
                    /* Copy the index in the columns */
                    osm->ir[toi] = tcia[ti].ir[nl] + cumsum(dc,i+1);
                    /* If necessary actualize the first element in the array */
                    if ( ts ){
                        osm->jc[nr+cumsum(dr,j+1)] = toi;
                        ts = 0;
                    }
                    toi ++;
                }
            }
            osm->jc[nr+cumsum(dr,j+1)+1] = toi;
        }
    }
    
}

int concatenate( const mxArray *ia, struct s_m *osm )
{
    mwSize *dc, *dr;
    mwSize i, nz, dn;
    const mwSize *ds;
    struct s_m *tcia;
    int cf = 0;

    ds = mxGetDimensions(ia);
    dn = mxGetNumberOfDimensions(ia);
    
    if ( dn > 2 ) {
        mexErrMsgTxt("Only two dimensional arrays are supported.");
    }

    if ( mxIsSparse( ia ) ){
        sparse2s_m(ia,osm);
        return 1;
    } else  if ( mxIsDouble(ia) ){
        full2sparse( ia, osm );
        return 1;
    } else if ( !mxIsCell(ia) ){
        mexErrMsgTxt("Input not supported.");
    } 
    
    /* If the input is a cell array */
    
    dc = mxCalloc(ds[0]+1,sizeof(mwSize));
    dr = mxCalloc(ds[1]+1,sizeof(mwSize));
    
    tcia = mxMalloc(sizeof(struct s_m) * ds[0] * ds[1]);

    nz = forwardPass(ia,tcia,dc,dr);

    /* Initilize the output */
    osm->nz = nz;
            
    osm->nc = cumsum(dc,ds[0]+1);
    osm->nr = cumsum(dr,ds[1]+1);

    osm->pr  = mxCalloc( nz, sizeof( double ) );
    osm->jc  = mxCalloc( (osm->nr)+1, sizeof( mwSize ) );
    osm->ir  = mxCalloc( nz, sizeof( mwSize ) );
    osm->jc[osm->nr] = nz;
    /* Check if there is any complex input */
    for ( i = 0; i < ds[0]*ds[1]; i++){
        if ( tcia[i].COMP ) {
            cf = 1;
            break;
        }
    }
    
    if ( cf ) {
        osm->pi  = mxCalloc(nz,sizeof(double));
        osm->COMP = 1;
    } else {
        osm->COMP = 0;        
    }

    /* Fill the output */

    backwardPass(tcia,osm,dc,dr,ds);

    /* Release memory */
    
    mxFree(dc);
    mxFree(dr);
    
    for ( i = 0; i < ds[0]*ds[1]; i++){
        free_s_m( tcia[i] );
    }
    mxFree( tcia );
    return 1;
     
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    if(nrhs < 1 || nrhs > 2) {
        mexErrMsgTxt("Incorrect number of inputs");
    } else if(nlhs>1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    if ( !( mxIsCell(prhs[0]) || mxIsNumeric(prhs[0]) || mxIsLogical(prhs[0]))){
        mexErrMsgTxt("Input is not supported.");
    }
    
    if ( nrhs == 1 ){
        double *pr, *pi = NULL;
        mwSize *ir, *jc;
        struct s_m osm;
        
        concatenate( prhs[0] , &osm );
        
        if ( osm.COMP ){
            plhs[0] = mxCreateSparse(osm.nc,osm.nr,osm.nz,mxCOMPLEX);
            pi      = mxGetPi(plhs[0]);
        } else {
            plhs[0] = mxCreateSparse(osm.nc,osm.nr,osm.nz,mxREAL);
        }
        
        pr = mxGetPr( plhs[0] );
        ir = mxGetIr( plhs[0] );
        jc = mxGetJc( plhs[0] );
        
        /* Transfer and clear memory */
        
        tmc_s_m(osm.COMP,osm.nr,osm.jc,osm.ir,osm.pr,osm.pi,
                jc,ir,pr,pi);
        free_s_m(osm);

        return;
        
    } else {
        mwSize i, j, dn, dd, td[2];
        const mwSize *ds;
        
        if ( !mxIsCell(prhs[0]) )
            mexErrMsgTxt("First argument should be a cell array");
        
        dn = mxGetNumberOfDimensions(prhs[0]);
        ds = mxGetDimensions(prhs[0]);
        
        if ( dn > 2 ) 
            mexErrMsgTxt("Only two dimensional cell arrays are supported.");
        
        if ( !mxIsNumeric(prhs[1]) )
            mexErrMsgTxt("Second argument should be numeric");
        
        dd = (mwSize ) mxGetScalar(prhs[1]);
        dd -= 1;
        
        if ( dd < 0 || dd > 1  )
            mexErrMsgTxt("Second argument should be equal to one or two.");
        
        
        /* Output cell array */
        
        td[!dd] = (mwSize ) ds[!dd];
        td[dd]  = (mwSize ) 1;
        
        plhs[0] = mxCreateCellArray(2,td);
                
        /* Iterate over dimensions */

        for ( i = 0 ; i < ds[!dd] ; i++){
            mwSize *jc, *ir, tds[2];
            double *pr, *pi = NULL;
            struct s_m osm;
            mxArray *ta;

            tds[dd] = ds[dd];
            tds[!dd] = 1;
      
            ta = mxCreateCellArray(2,tds);
      
            
            for ( j = 0 ; j < ds[dd] ; j++){
                mxArray *tsa;
                mwSize ii[2];
                
                ii[0] = j;
                ii[1] = i;
                
                tsa = mxGetCell(prhs[0],ii[dd]+(ii[!dd]*ds[0]));
                
                if (tsa != NULL)
                    mxSetCell(ta,j,mxDuplicateArray(tsa));
            }

            concatenate( ta, &osm );

            if ( osm.COMP ) {
                mxSetCell(plhs[0],i,mxCreateSparse(osm.nc,osm.nr,osm.nz,mxCOMPLEX));
                pi = mxGetPi( mxGetCell(plhs[0],i) );
            } else {
                mxSetCell(plhs[0],i, mxCreateSparse(osm.nc,osm.nr,osm.nz,mxREAL) );
            }
            
            pr = mxGetPr( mxGetCell(plhs[0],i) );
            ir = mxGetIr( mxGetCell(plhs[0],i) );
            jc = mxGetJc( mxGetCell(plhs[0],i) );
                       
            /* Transfer and clear memory */

            tmc_s_m(osm.COMP,osm.nr,osm.jc,osm.ir,osm.pr,osm.pi,jc,ir,pr,pi);
            free_s_m(osm);
            mxDestroyArray(ta);
        }
    return;
    }
}
 
