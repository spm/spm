#include "mex.h"
#include "shoot_openmp.h"
#ifdef _OPENMP
#   include <omp.h>
#endif

static int num_threads = 1;

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
        num_threads = 1;
#   endif
}

int get_num_threads()
{
    return(num_threads);
}

int get_num_procs()
{
#ifdef _OPENMP
    return(omp_get_num_procs());
#else
    return(1);
#endif
}

