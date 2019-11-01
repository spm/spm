#include "shoot_openmp.h"
#ifdef _OPENMP
#   include <omp.h>
#endif
#ifndef SPM_NUMTHREADS_DEFAULT
#   ifdef _OPENMP
#       define SPM_NUMTHREADS_DEFAULT 1
#   else
#       define SPM_NUMTHREADS_DEFAULT 0
#   endif
#endif

static int num_threads = SPM_NUMTHREADS_DEFAULT;

void set_num_threads(int t)
{
#   ifdef _OPENMP
        if(t == 0)
            num_threads = 1;
        else if(t > 0)
            num_threads = t;
        else
            num_threads = -t * get_num_procs();
        omp_set_num_threads(num_threads);
#   else
        num_threads = 0;
#   endif
}

int get_num_threads()
{
#ifdef _OPENMP
    return(num_threads);
#else
    return(0);
#endif
}

int get_num_procs()
{
#ifdef _OPENMP
    return(omp_get_num_procs());
#else
    return(0);
#endif
}

