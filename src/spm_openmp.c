/*
 * $Id: spm_openmp.c 7689 2019-11-07 13:25:47Z guillaume $
 */

#include "spm_openmp.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef SPM_WIN32
#include <windows.h>
#endif

static char SPM_NUM_THREADS[] = "SPM_NUM_THREADS";

void spm_set_num_threads(int t)
{
    int num_threads;
    char spm_num_threads[8];
#ifdef _OPENMP
    if(t == 0)
        num_threads = 1;
    else if(t > 0)
        num_threads = t;
    else
        num_threads = -t * omp_get_num_procs();
    omp_set_num_threads(num_threads);
#else
    num_threads = 1;
#endif
    sprintf(spm_num_threads,"%d",num_threads);
#ifdef SPM_WIN32
    _putenv_s(SPM_NUM_THREADS,spm_num_threads);
#else
    setenv(SPM_NUM_THREADS,spm_num_threads,1);
#endif
}

int spm_get_num_threads()
{
#ifdef _OPENMP
    
#ifdef SPM_WIN32
    static char num_threads[MAX_PATH];
    const DWORD ret = GetEnvironmentVariable(SPM_NUM_THREADS,num_threads,MAX_PATH);
    int sts = (ret == 0) | (ret > MAX_PATH);
#else
    char *num_threads = getenv(SPM_NUM_THREADS);
    int sts = num_threads == NULL;
#endif
    if (sts)
        return(1); /* default: alternative is -1 */
    else
        return(atoi(num_threads));
    
#else
    return(1);
#endif
}
