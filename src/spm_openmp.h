/*
 * $Id: spm_openmp.h 7688 2019-11-07 12:24:42Z guillaume $
 */

#if ! defined (__SPM_OPENMP_H_)
#define __SPM_OPENMP_H_

#ifdef _OPENMP
#include <omp.h>
#define SPM_OMP_PRAGMA(x) _Pragma (#x)
#else
#define SPM_OMP_PRAGMA(x)
#endif

extern void spm_set_num_threads(int t);
extern int  spm_get_num_threads();

#endif
