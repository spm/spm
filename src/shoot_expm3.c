/* $Id: shoot_expm3.c 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */

#include <mex.h>
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
static float *sub33(float *a, float *b, float *c)
{
    int i;
    for(i=0; i<9; i++)
        c[i] = a[i] - b[i];
    return(c);
}

static float *add33(float *a, float *b, float *c)
{
    int i;
    for(i=0; i<9; i++)
        c[i] = a[i] + b[i];
    return(c);
}

static float *scale33(float *a, float s, float *b)
{
    int i;
    for(i=0; i<9; i++)
        b[i] = a[i]*s;
    return(b);
}

static float *eye33(float *a)
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
    return(a);
}

static float *mul33(float *a, float *b, float *c)
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
    return(c);
}

static float *div33(float *a, float *b, float *c)
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
    return(c);
}

static void pade33(float *a, float *l)
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

static void pade22(float *a, float *l)
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
    rm = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);
    r  = fabs(a[3]) + fabs(a[4]) + fabs(a[5]);
    if (r>rm) rm = r;
    r  = fabs(a[6]) + fabs(a[7]) + fabs(a[8]);
    if (r>rm) rm = r;
    return(rm);
}

static float *assign33(float *a, float *b)
{
    int i;
    for(i=0; i<9; i++)
        b[i] = a[i];
    return(b);
}

void expm33(float *a, float *l)
{
    /* See expm.m in MATLAB or http://mathworld.wolfram.com/PadeApproximant.html */
    int K;
    K = (int)ceil(log((double)(norm1(a)*2.3481))*1.44269504088896);
    if (K>0)
    {
        float b[9];
        float s = 1.0/pow2(K);
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

void expm22(float *a, float *l)
{
    /* See expm.m in MATLAB or http://mathworld.wolfram.com/PadeApproximant.html */
    int K;
    K = (int)ceil(log((double)(norm1(a)*12.356))*1.44269504088896);
    if (K>0)
    {
        float b[9];
        float s = 1.0/pow2(K);
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


