/*
 * $Id: spm_mex.h 7686 2019-11-06 13:18:06Z guillaume $
 */

#if defined (MATLAB_MEX_FILE) || defined (HAVE_OCTAVE)
#   include "mex.h"
#else
#   if ! defined (__SPM_MEX_H_)
#       define __SPM_MEX_H_
#       include <stddef.h>
#       include <stdio.h>
#       include <math.h>
#       define  mwSignedIndex signed   long long int
#       define  mwIndex       unsigned long long int
#       define  mwSize        unsigned long long int
#       define  mxIsFinite(x) isfinite(x)
#       define  mxGetNaN() (NAN)
#       define  mexErrMsgTxt(msg) perror(msg)
#   endif
#endif
