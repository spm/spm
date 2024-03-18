/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include "cuheader.h"

/*
static float wt0(float x)
{
    return((fabsf(x) > 0.5f) ? 0.0f : 1.0f);
}

static float wt1(float x)
{
    x = fabsf(x);
    return((x > 1.0f) ? 0.0f : (1.0f - x));
}
*/

__device__ float wt2(float x)
{
    x = fabsf(x);
    if(x < 0.5f)
        return(0.75f - x*x);
    if(x < 1.5f)
    {
        x = 1.5f - x;
        return(0.5f*x*x);
    }
    return(0.0f);
}

__device__ float wt3(float x)
{
    x = fabsf(x);
    if(x < 1.0f)
        return(x*x*(x - 2.0f)*0.5f + 2.0f/3.0f);
    if(x < 2.0f)
    {
        x = 2.0f - x;
        return(x*x*x*(1.0f/6.0f));
    }
    return(0.0f);
}

__device__ float wt4(float x)
{
    x = fabsf(x);
    if(x < 0.5f)
    {
        x *= x;
        return(x*(x*0.25f - 0.625f) + 115.0f/192.0f);
    }
    if(x < 1.5f)
        return(x*(x*(x*(5.0f/6.0f - x*(1.0f/6.0f)) - 1.25f) + 5.0f/24.0f) + 55.0f/96.0f);
    if(x < 2.5f)
    {
        x -= 2.5f;
        x *= x;
        return(x*x*(1.0f/24.0f));
    }
    return(0.0);
}


/* Note that dp[i] is 1 more than the interpolation degree */
__device__ SSIZE_t weights(const USIZE_t d, float x, /*@OUT@*/float w[])
{
    USIZE_t k;
    SSIZE_t i = (SSIZE_t)ceilf(x-0.5f*(float)d);
    x -= (float)i;

    switch (d){
    case 2:
        w[0] = 1.0f-x;
        w[1] = x;
        break;
    case 1:
        w[0] = 1.0f;
        break;
    case 3:
        for(k=0; k<=2; k++) w[k] = wt2(x-(float)k);
        break;
    case 4:
        for(k=0; k<=3; k++) w[k] = wt3(x-(float)k);
        break;
    case 5:
        for(k=0; k<=4; k++) w[k] = wt4(x-(float)k);
        break;
    default:
        w[0] = 1.0f-x;
        w[1] = x;
    }
    return i;
}


__device__ float dwt2(float x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabsf(x);

    if(x < 0.5f)
        return(-2*x*s);
    if(x < 1.5f)
        return((x - 1.5f)*s);
    return(0.0f);
}

__device__ float dwt3(float x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabsf(x);

    if(x < 1.0f)
        return(x*(1.5f*x - 2.0f)*s);
    if(x < 2.0f)
    {
        x = x - 2.0f;
        return(-0.5f*x*x*s);
    }
    return(0.0f);
}

__device__ float dwt4(float x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabsf(x);

    if(x < 0.5f)
        return((x*(x*x - 5.0f/4.0f))*s);
    if(x < 1.5f)
        return((x*(x*(x*(-2.0f/3.0f) + 2.5f) - 5.0f/2.0f) + 5.0f/24.0f)*s);
    if(x < 2.5f)
    {
        x = x*2.0f - 5.0f;
        return((1.0f/48.0f)*x*x*x*s);
    }
    return(0.0f);
}


/* Note that dp[i] is 1 more than the interpolation degree */
__device__ SSIZE_t dweights(const USIZE_t d, float x, /*@OUT@*/float w[])
{
    USIZE_t k;
    SSIZE_t i = (SSIZE_t)ceilf(x-0.5f*(float)d);
    x -= (float)i;

    switch (d){
    case 2:
        w[0] = -1.0f;
        w[1] =  1.0f;
        break;
    case 1:
        w[0] = 0.0f;
        break;
    case 3:
        for(k=0; k<=2; k++) w[k] = dwt2(x-(float)k);
        break;
    case 4:
        for(k=0; k<=3; k++) w[k] = dwt3(x-(float)k);
        break;
    case 5:
        for(k=0; k<=4; k++) w[k] = dwt4(x-(float)k);
        break;
    default:
        w[0] = -1.0f;
        w[1] =  1.0f;
    }
    return i;
}

/*
__device__ float hwt2(float x)
{
    x = fabsf(x);

    if(x < 0.5f)
        return(-2.0f);
    if(x < 1.5f)
        return(1.0f);
    return(0.0f);
}
*/

__device__ float hwt3(float x)
{
    x = fabsf(x);

    if(x < 1.0f)
        return(3.0f*x - 2.0f);
    if(x < 2.0f)
        return(2.0f - x);
    return(0.0f);
}

__device__ float hwt4(float x)
{
    x = fabsf(x);

    if(x < 0.5f)
        return(3.0f*x*x - 5.0f/4.0f);
    if(x < 1.5f)
        return(x*(-2.0f*x + 5.0f) - 5.0f/2.0f);
    if(x < 2.5f)
    {
        x = x*2.0f - 5.0f;
        return(x*x/8.0f);
    }
    return(0.0f);
}

__device__ SSIZE_t hweights(const USIZE_t d, float x, /*@OUT@*/float w[])
{
    USIZE_t k;
    SSIZE_t i = (SSIZE_t)ceilf(x-0.5f*(float)d);
    x -= (float)i;

    switch (d){
    case 2:
        w[0] = 0.0f;
        w[1] = 0.0f;
        break;
    case 1:
        w[0] = 0.0f;
        break;
    case 3:
     /* for(k=0; k<=2; k++) w[k] = hwt2(x-(float)k); */
        w[0] =  1.0f;
        w[1] = -2.0f;
        w[2] =  1.0f;
        break;
    case 4:
        for(k=0; k<=3; k++) w[k] = hwt3(x-(float)k);
        break;
    case 5:
        for(k=0; k<=4; k++) w[k] = hwt4(x-(float)k);
        break;
    default:
        w[0] = 0.0f;
        w[1] = 0.0f;
    }
    return i;
}
