/*
 * $Id: spm_existfile.c 6428 2015-05-06 14:09:04Z guillaume $
 * Guillaume Flandin
 */

#ifndef MATLAB_MEX_FILE
#undef  _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#include <sys/stat.h>
#if defined(_WIN32)
#define structStat struct _stati64
#define getFileStat _stati64
#elif defined(__APPLE__)
#define structStat struct stat
#define getFileStat stat
#else
#define structStat struct stat64
#define getFileStat stat64
#endif
#else
#include "io64.h"
#endif
#include "mex.h"

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int status     = 0;
    char *filename = NULL;
    structStat stbuf;
    
    if (nrhs != 1)
    {
        mexErrMsgTxt("One input only required.");
    }
    else
    {
        if (!mxIsChar(prhs[0]))
        {
            mexErrMsgTxt("Input must be a string.");
        }
        filename = mxArrayToString(prhs[0]);
        
        if ((getFileStat(filename, &stbuf) == 0) && (S_ISREG(stbuf.st_mode)))
        {
            status = 1;
        }

        mxFree(filename);
    }
        
    plhs[0] = mxCreateLogicalScalar(status);
}
