#ifdef MATLAB_MEX_FILE
#   include "mex.h"
#else
#   include <stddef.h>
#   include <math.h>
#   define  mxIsfinite(x) isfinite(x)
#   define  mxGetNaN() (NAN)
#   define  mwSignedIndex signed   long long int
#   define  mwIndex       unsigned long long int
#   define  mwSize        unsigned long long int
#endif

