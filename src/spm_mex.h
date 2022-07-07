/*
 * John Ashburner
 * Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging
 */

#if defined (MATLAB_MEX_FILE) || defined (OCTAVE_MEX_FILE)
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
#       define  mxIsFinite(x) (bool)isfinite(x)
#       define  mxGetNaN() (NAN)
#       define  mexErrMsgTxt(msg) perror(msg)
#   endif
#endif
