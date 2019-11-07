/*
 * $Id: spm_openmp.c 7687 2019-11-07 11:26:02Z guillaume $
 */

#include "spm_openmp.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
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
    setenv(SPM_NUM_THREADS,spm_num_threads,1);
}

int spm_get_num_threads()
{
#ifdef _OPENMP
    char *num_threads = getenv(SPM_NUM_THREADS);
    if (num_threads == NULL)
        return(1); /* default: alternative is -1 */
    else
        return(atoi(num_threads));
#else
    return(1);
#endif
}
