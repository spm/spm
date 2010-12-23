/*
 * $Id: spm_existfile.c 4145 2010-12-23 15:18:30Z guillaume $
 * Guillaume Flandin
 */
 
#define _FILE_OFFSET_BITS 64

#include "mex.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int status     = 0;
    char *filename = NULL;
    FILE *fid      = NULL;
    
    if (nrhs != 1)
        mexErrMsgTxt("One input only required.");
    else
    {
        if (!mxIsChar(prhs[0]))
            mexErrMsgTxt("Input must be a string.");
        filename = mxArrayToString(prhs[0]);
        fid = fopen(filename,"r");
        if (fid != NULL)
        {
            status = 1;
            fclose(fid);
        }
        mxFree(filename);
    }
        
    plhs[0] = mxCreateLogicalScalar(status);
}
