/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#ifndef _CUHDR
#    define _CUHDR
#    include<stddef.h>
#    ifndef CUDA
#        define __device__ static
#        define __global__
#        include<math.h>
#        define atomicAdd(p,v) (*p)+=(v)
#    endif
#    define USIZE_t size_t
#    define SSIZE_t ptrdiff_t
/*
#    define USIZE_t unsigned int
#    define SSIZE_t int
*/
#endif
