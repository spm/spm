/* $Id: shoot_expm3.c 7700 2019-11-21 17:09:15Z john $ */
/* (c) John Ashburner (2011) */

#include "spm_mex.h"
#include <math.h>
#include <stdio.h>

extern double   log(double x);
extern double   exp(double x);

static int pow2(int k)
{
    int j0, td = 1;
    for(j0=0; j0<k; j0++)
        td = td*2;
    return(td);
}

/* syms a[0] a[1] a[2] a[3] a[4] a[5] a[6] a[7] a[8]
   syms b[0] b[1] b[2] b[3] b[4] b[5] b[6] b[7] b[8]
   A = [a[0] a[1] a[2]; a[3] a[4] a[5]; a[6] a[7] a[8]].';
   B = [b[0] b[1] b[2]; b[3] b[4] b[5]; b[6] b[7] b[8]].';
   C = A*B;
   C = C(:)

   C = A\B;
   C = C(:)
*/
static void sub33(float *a, float *b, /*@out@*/float *c)
{
    int i;
    for(i=0; i<9; i++)
        c[i] = a[i] - b[i];
}

static void add33(float *a, float *b, /*@out@*/float *c)
{
    int i;
    for(i=0; i<9; i++)
        c[i] = a[i] + b[i];
}

static void scale33(float *a, float s, /*@out@*/float *b)
{
    int i;
    for(i=0; i<9; i++)
        b[i] = a[i]*s;
}

static void eye33(/*@out@*/float *a)
{
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 0.0;
    a[3] = 0.0;
    a[4] = 1.0;
    a[5] = 0.0;
    a[6] = 0.0;
    a[7] = 0.0;
    a[8] = 1.0;
}

static void mul33(float *a, float *b, /*@out@*/float *c)
{
    c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
    c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
    c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
    c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
    c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
    c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
    c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
    c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
    c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
}

static void div33(float *a, float *b, /*@out@*/float *c)
{
    float d = a[0]*(a[4]*a[8] - a[5]*a[7]) + a[1]*(a[5]*a[6] - a[3]*a[8]) + a[2]*(a[3]*a[7] - a[4]*a[6]);
    c[0] =   (a[3]*(a[7]*b[2] - a[8]*b[1]) + a[4]*(a[8]*b[0] - a[6]*b[2]) + a[5]*(a[6]*b[1] - a[7]*b[0]))/d;
    c[1] =  -(a[0]*(a[7]*b[2] - a[8]*b[1]) + a[1]*(a[8]*b[0] - a[6]*b[2]) + a[2]*(a[6]*b[1] - a[7]*b[0]))/d;
    c[2] =   (a[0]*(a[4]*b[2] - a[5]*b[1]) + a[1]*(a[5]*b[0] - a[3]*b[2]) + a[2]*(a[3]*b[1] - a[4]*b[0]))/d;
    c[3] =   (a[3]*(a[7]*b[5] - a[8]*b[4]) + a[4]*(a[8]*b[3] - a[6]*b[5]) + a[5]*(a[6]*b[4] - a[7]*b[3]))/d;
    c[4] =  -(a[0]*(a[7]*b[5] - a[8]*b[4]) + a[1]*(a[8]*b[3] - a[6]*b[5]) + a[2]*(a[6]*b[4] - a[7]*b[3]))/d;
    c[5] =   (a[0]*(a[4]*b[5] - a[5]*b[4]) + a[1]*(a[5]*b[3] - a[3]*b[5]) + a[2]*(a[3]*b[4] - a[4]*b[3]))/d;
    c[6] =   (a[3]*(a[7]*b[8] - a[8]*b[7]) + a[4]*(a[8]*b[6] - a[6]*b[8]) + a[5]*(a[6]*b[7] - a[7]*b[6]))/d;
    c[7] =  -(a[0]*(a[7]*b[8] - a[8]*b[7]) + a[1]*(a[8]*b[6] - a[6]*b[8]) + a[2]*(a[6]*b[7] - a[7]*b[6]))/d;
    c[8] =   (a[0]*(a[4]*b[8] - a[5]*b[7]) + a[1]*(a[5]*b[6] - a[3]*b[8]) + a[2]*(a[3]*b[7] - a[4]*b[6]))/d;
}

static void pade33(float *a, /*@out@*/ float *l)
{
    float u[9], v[9], num[9], den[9], a0[9], a2[9], a3[9];
    
    eye33(a0);
    mul33(a,a,a2);
    mul33(a2,a,a3);
    scale33(a0,120.0,a0);
    scale33(a2, 12.0,a2);
    add33(a0,a2,u);
    scale33(a,60.0,v);
    add33(v,a3,v);
    
    add33(u,v,num);
    sub33(u,v,den);
    div33(den,num,l);
}

static void pade22(float *a, /*@out@*/ float *l)
{
    float u[9], v[9], num[9], den[9], a0[9], a2[9];
    
    eye33(a0);
    mul33(a,a,a2);
    scale33(a0,12.0,a0);
    add33(a0,a2,u);
    scale33(a,6.0,v);
    
    add33(u,v,num);
    sub33(u,v,den);
    div33(den,num,l);
}

static float norm1(float *a)
{
    float r, rm;
    rm = (float)fabs(a[0]) + (float)fabs(a[1]) + (float)fabs(a[2]);
    r  = (float)fabs(a[3]) + (float)fabs(a[4]) + (float)fabs(a[5]);
    if (r>rm) rm = r;
    r  = (float)fabs(a[6]) + (float)fabs(a[7]) + (float)fabs(a[8]);
    if (r>rm) rm = r;
    return(rm);
}

static  void assign33(float *a, float *b)
{
    int i;
    for(i=0; i<9; i++)
        b[i] = a[i];
}

void expm33(float *a, /*@out@*/ float *l)
{
    /* See expm.m in MATLAB or http://mathworld.wolfram.com/PadeApproximant.html */
    int K;
    K = (int)ceil(log((double)(norm1(a)*2.3481))*1.44269504088896);
    if (K>0)
    {
        float b[9];
        float s = 1.0f/(float)pow2(K);
        int i;
        scale33(a,s,b);
        pade33(b, l);
        for(i=0; i<K; i++)
        {
            assign33(l,b);
            mul33(b,b,l);
        }
    }
    else
        pade33(a, l);
}

void expm22(float *a, /*@out@*/ float *l)
{
    /* See expm.m in MATLAB or http://mathworld.wolfram.com/PadeApproximant.html */
    int K;
    K = (int)ceil(log((double)(norm1(a)*12.356))*1.44269504088896);
    if (K>0)
    {
        float b[9];
        float s = 1.0f/(float)pow2(K);
        int i;
        scale33(a,s,b);
        pade22(b, l);
        for(i=0; i<K; i++)
        {
            assign33(l,b);
            mul33(b,b,l);
        }
    }
    else
        pade22(a, l);
}


