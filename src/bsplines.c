/*
 * $Id: bsplines.c 4624 2012-01-13 13:27:08Z john $
 * John Ashburner
 */
 
/*
 * This code is a modified version of that of Philippe Thevenaz, which I took from:
 *  http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs. Many thanks to Philippe Thevenaz for advice with the code.
 *
 * See:
 *  M. Unser, A. Aldroubi and M. Eden.
 *  "B-Spline Signal Processing: Part I-Theory,"
 *  IEEE Transactions on Signal Processing 41(2):821-832 (1993).
 *
 *  M. Unser, A. Aldroubi and M. Eden.
 *  "B-Spline Signal Processing: Part II-Efficient Design and Applications,"
 *  IEEE Transactions on Signal Processing 41(2):834-848 (1993).
 *
 *  M. Unser.
 *  "Splines: A Perfect Fit for Signal and Image Processing,"
 *  IEEE Signal Processing Magazine 16(6):22-38 (1999).
 *
*/
#include "mex.h"
#include <stdlib.h>
#include <math.h>
#ifdef IMAGE_SINGLE
#define IMAGE_DTYPE float
#else
#define IMAGE_DTYPE double
#endif

/***************************************************************************************
Starting periodic boundary condition based on Eq. 2.6 of Unser's 2nd 1993 paper.
    c - vector of unfiltered data
    m - length of c
    p - pole (root of polynomial)
    function returns value that c[0] should initially take

The expression for the first pass of the recursive convolution is:
    for (i=1; i<m; i++) c[i] += p*c[i-1];

If m==4, then:
    c0 = c0 + p*c3;
    c1 = c1 + p*c0;
    c2 = c2 + p*c1
    c3 = c3 + p*c2;
    c0 = c0 + p*c3;
    etc...
After recursive substitution, c0 becomes:
    (1  +p^4+p^8 +p^12 ...)*c0 +
    (p  +p^5+p^9 +p^13 ...)*c3 +
    (p^2+p^6+p^10+p^14 ...)*c2 +
    (p^3+p^7+p^11+p^15 ...)*c1

Using maple...
    sum('p^(k*m+n)','k'=0..infinity)
These series converge to...
    (p^n)/(1-p^m)

So c0 becomes:
    (c0 + c3*p + c2*p^2 + c1*p^3)/(1-p^4)
*/
static double cc_wrap(IMAGE_DTYPE c[], int m, double p)
{
    double s, pi;
    int    i, m1;

    m1 = ceil(-30/log(fabs(p)));
    if (m1>m) m1=m;

    pi   = p;
    s    = c[0];
    for (i=1; i<m1; i++)
    {
        s   += pi*c[m-i];
        pi  *= p;
    }
    return(s/(1.0-pi));
}

/***************************************************************************************
Starting mirrored boundary condition based on Eq. 2.6 of Unser's 2nd 1993 paper.
    c - vector of unfiltered data
    m - length of c
    p - pole (root of polynomial)
    function returns value that c[0] should initially take
*/
static double cc_mirror(IMAGE_DTYPE c[], int m, double p)
{
    double s, pi, p2i, ip;
    int    i, m1;

    m1 = ceil(-30/log(fabs(p)));
    if (m1 < m)
    {
        pi = p;
        s  = c[0];
        for (i=1; i<m1; i++)
        {
            s  += pi * c[i];
            pi *= p;
        }
        return(s);
    }
    else
    {
        pi   = p;
        ip   = 1.0/p;
        p2i  = pow(p,m-1.0);
        s    = c[0] + p2i*c[m-1];
        p2i *= p2i * ip;
        for (i=1; i<m-1; i++)
        {
            s   += (pi+p2i)*c[i];
            pi  *= p;
            p2i *= ip;
        }
        return(s/(1.0-pi*pi));
    }
}

/***************************************************************************************
End periodic boundary condition
    c - first pass filtered data
    m - length of filtered data (must be > 1)
    p - pole
    function returns value for c[m-1] before 2nd filter pass

The expression for the second pass of the recursive convolution is:
    for (i=m-2; i>=0; i--) c[i] = p*(c[i+1]-c[i]);
If m==4, then:
    c3 = p*(c0-c3);
    c2 = p*(c3-c2);
    c1 = p*(c2-c1);
    c0 = p*(c1-c0);
    c3 = p*(c0-c3);
    etc...

After recursive substitution, c0 becomes:
    -(p  +p^5+p^9  ...)*c3
    -(p^2+p^6+p^10 ...)*c0
    -(p^3+p^7+p^11 ...)*c1
    -(p^4+p^8+p^12 ...)*c2

These series converge to...
    (p^n)/(p^m-1)

So c0 becomes:
    (c3*p + c0*p^2 + c1*p^3 + c2*p^4)/(p^4-1)
*/
static double icc_wrap(IMAGE_DTYPE c[], int m, double p)
{
    double s, pi;
    int    i, m1;

    m1 = ceil(-30/log(fabs(p)));
    if (m1>m) m1=m;

    pi = p;
    s  = pi*c[m-1];
    for (i=0; i<m1-1; i++)
    {
        pi  *= p;
        s   += pi*c[i];
    }
    return(s/(pi-1.0));
}

/***************************************************************************************
End mirrored boundary condition
    c - first pass filtered data
    m - length of filtered data (must be > 1)
    p - pole
    function returns value for c[m-1] before 2nd filter pass
*/
static double icc_mirror(IMAGE_DTYPE c[],int m, double p)
{
    return((p/(p*p-1.0))*(p*c[m-2]+c[m-1]));
}

/***************************************************************************************
Compute gains required for zero-pole representation - see tf2zp.m in Matlab's
 Signal Processing Toolbox.
    p - poles
    np - number of poles
    function returns the gain of the system
*/
static double gain(double p[], int np)
{
    int j;
    double lambda = 1.0;
    for (j = 0; j < np; j++)
        lambda = lambda*(1.0-p[j])*(1.0-1.0/p[j]);
    return(lambda);
}

/***************************************************************************************
One dimensional recursive filtering - assuming wrapped boundaries
See Eq. 2.5 of Unsers 2nd 1993 paper.
    c - original vector on input, coefficients on output
    m - length of vector
    p - poles (polynomial roots)
    np - number of poles
*/
void splinc_wrap(IMAGE_DTYPE c[], int m, double p[], int np)
{
    double lambda = 1.0;
    int i, k;

    if (m == 1) return;

    /* compute gain and apply it */
    lambda = gain(p,np);
    for (i = 0; i < m; i++)
        c[i] *= lambda;

    /* loop over poles */
    for (k = 0; k < np; k++)
    {
        double pp = p[k];
        c[0] = cc_wrap(c, m, pp);

        for (i=1; i<m; i++)
            c[i] += pp*c[i-1];

        c[m-1] = icc_wrap(c, m, pp);
        for (i=m-2; i>=0; i--)
            c[i] = pp*(c[i+1]-c[i]);
    }
}

/***************************************************************************************
One dimensional recursive filtering - assuming mirror boundaries
See Eq. 2.5 of Unsers 2nd 1993 paper.
    c - original vector on input, coefficients on output
    m - length of vector
    p - poles (polynomial roots)
    np - number of poles
*/
void splinc_mirror(IMAGE_DTYPE c[], int m, double p[], int np)
{
    double lambda = 1.0;
    int i, k;

    if (m == 1) return;

    /* compute gain and apply it */
    lambda = gain(p,np);
    for (i = 0; i < m; i++)
        c[i] *= lambda;

    /* loop over poles */
    for (k = 0; k < np; k++)
    {
        double pp = p[k];
        c[0] = cc_mirror(c, m, pp);

        for (i=1; i<m; i++)
            c[i] += pp*c[i-1];

        c[m-1] = icc_mirror(c, m, pp);
        for (i=m-2; i>=0; i--)
            c[i] = pp*(c[i+1]-c[i]);
    }
}

/***************************************************************************************
Return roots of B-spline kernels.
     d - degree of B-spline
     np - number of roots of magnitude less than one
     p - roots.
*/
int get_poles(int d, int *np, double p[])
{
    /* Return polynomial roots that are less than one. */
    switch (d) {
        case 0:
            *np = 0;
            break;
        case 1:
            *np = 0;
            break;
        case 2: /* roots([1 6 1]) */
            *np = 1;
            p[0] = sqrt(8.0) - 3.0;
            break;
        case 3: /* roots([1 4 1]) */
            *np = 1;
            p[0] = sqrt(3.0) - 2.0;
            break;
        case 4: /* roots([1 76 230 76 1]) */
            *np = 2;
            p[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
            p[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
            break;
        case 5: /* roots([1 26 66 26 1]) */
            *np   = 2;
            p[0] = sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25) - 6.5;
            p[1] = sqrt(67.5 + sqrt(4436.25)) - sqrt(26.25) - 6.5;
            break;
        case 6: /* roots([1 722 10543 23548 10543 722 1]) */
            *np   = 3;
            p[0] = -0.488294589303044755130118038883789062112279161239377608394;
            p[1] = -0.081679271076237512597937765737059080653379610398148178525368;
            p[2] = -0.00141415180832581775108724397655859252786416905534669851652709;
            break;
        case 7: /* roots([1 120 1191 2416 1191 120 1]) */
            *np   = 3;
            p[0] = -0.5352804307964381655424037816816460718339231523426924148812;
            p[1] = -0.122554615192326690515272264359357343605486549427295558490763;
            p[2] = -0.0091486948096082769285930216516478534156925639545994482648003;
            break;
        default:
            return(1);
    }
    return(0);
}
/***************************************************************************************/

/***************************************************************************************
Different degrees of B-splines
    x - position relative to origin
    returns value of basis function at x
*/

/*static double wt1(double x)
{
    x = fabs(x);
    return((x > 1.0) ? (0.0) : (1.0 - x));
}*/

static double wt2(double x)
{
    x = fabs(x);
    if (x < 0.5)
        return(0.75 - x*x);
    if (x < 1.5)
    {
        x = 1.5 - x;
        return(0.5*x*x);
    }
    return(0.0);
}

static double wt3(double x)
{
    x = fabs(x);
    if (x < 1.0)
        return(x*x*(x - 2.0)*0.5 + 2.0/3.0);
    if (x < 2.0)
    {
        x = 2.0 - x;
        return(x*x*x*(1.0/6.0));
    }
    return(0.0);
}

static double wt4(double x)
{
    x = fabs(x);
    if (x < 0.5)
    {
        x *= x;
        return(x*(x*0.25 - 0.625) + 115.0/192.0);
    }
    if (x < 1.5)
        return(x*(x*(x*(5.0/6.0 - x*(1.0/6.0)) - 1.25) + 5.0/24.0) + 55.0/96.0);
    if (x < 2.5)
    {
        x -= 2.5;
        x *= x;
        return(x*x*(1.0/24.0));
    }
    return(0.0);
}

static double wt5(double x)
{
    double y;
    x = fabs(x);
    if (x < 1.0)
    {
        y = x*x;
        return(y*(y*(0.25 - x*(1.0/12.0)) - 0.5) + 0.55);
    }
    if (x < 2.0)
        return(x*(x*(x*(x*(x*(1.0/24.0) - 0.375) + 1.25) - 1.75) + 0.625) + 0.425);
    if (x < 3.0)
    {
        y = 3.0 - x;
        x = y*y;
        return(y*x*x*(1.0/120.0));
    }
    return(0.0);
}

static double wt6(double x)
{
    x = fabs(x);
    if (x < 0.5)
    {
        x *= x;
        return(x*(x*(7.0/48.0 - x*(1.0/36.0)) - 77.0/192.0) + 5887.0/11520.0);
    }
    if (x < 1.5)
        return(x*(x*(x*(x*(x*(x*(1.0/48.0) - 7.0/48.0) + 0.328125)
             - 35.0/288.0) - 91.0/256.0) - 7.0/768.0) + 7861.0/15360.0);
    if (x < 2.5)
        return(x*(x*(x*(x*(x*(7.0/60.0 - x*(1.0/120.0)) - 0.65625)
            + 133.0/72.0) - 2.5703125) + 1267.0/960.0) + 1379.0/7680.0);
    if (x < 3.5)
    {
        x -= 3.5;
        x *= x*x;
        return(x*x*(1.0/720.0));
    }
    return(0.0);
}

static double wt7(double x)
{
    double y;

    x = fabs(x);
    if (x < 1.0)
    {
        y = x*x;
        return(y*(y*(y*(x*(1.0/144.0) - 1.0/36.0) + 1.0/9.0) - 1.0/3.0)
            + 151.0/315.0);
    }
    if (x < 2.0)
        return(x*(x*(x*(x*(x*(x*(0.05 - x*(1.0/240.0)) - 7.0/30.0) + 0.5)
            - 7.0/18.0) - 0.1) - 7.0/90.0) + 103.0/210.0);
    if (x < 3.0)
        return(x*(x*(x*(x*(x*(x*(x*(1.0/720.0) - 1.0/36.0) + 7.0/30.0)
            - 19.0/18.0) + 49.0/18.0) - 23.0/6.0) + 217.0/90.0) - 139.0/630.0);
    if (x < 4.0)
    {
        y = 4.0 - x;
        x = y*y*y;
        return(x*x*y*(1.0/5040.0));
    }
    return(0.0);
}

/***************************************************************************************
Derivatives of different degrees of B-splines
    x - position relative to origin
    returns derivative of basis function at x
*/

static double dwt2(double x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabs(x);

    if (x < 0.5)
        return(-2*x*s);
    if (x < 1.5)
        return((x - 1.5)*s);
    return(0.0);
}

static double dwt3(double x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabs(x);


    if (x < 1.0)
        return(x*(1.5*x - 2.0)*s);
    if (x < 2.0)
    {
        x = x - 2.0;
        return(-0.5*x*x*s);
    }
    return(0.0);
}

static double dwt4(double x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabs(x);

    if (x < 0.5)
    {
        return((x*(x*x - 5.0/4.0))*s);
    }
    if (x < 1.5)
        return((x*(x*(x*(-2.0/3.0) + 2.5) - 5.0/2.0) + 5.0/24.0)*s);
    if (x < 2.5)
    {
        x = x*2.0 - 5.0;
        return((1.0/48.0)*x*x*x*s);
    }
    return(0.0);
}

static double dwt5(double x)
{
    int s;
    s = (x>0 ? 1 : -1);
    x = fabs(x);

    if (x < 1.0)
        return((x*(x*(x*(x*(-5.0/12.0) + 1.0)) - 1.0))*s);
    if (x < 2.0)
        return((x*(x*(x*(x*(5.0/24.0) - 1.5) + 3.75) - 3.5) + 0.625)*s);
    if (x < 3.0)
    {
        x -= 3.0;
        x *= x;
        return((-1.0/24.0)*x*x*s);
    }
    return(0.0);
}

static double dwt6(double x)
{
    double y;
    int s;
    s = (x>0 ? 1 : -1);
    x = fabs(x);

    if (x < 0.5)
    {
        y = x*x;
        return(x*((7.0/12.0)*y - (1.0/6.0)*y*y - (77.0/96.0))*s);
    }
    if (x < 1.5)
        return((x*(x*(x*(x*(x*0.125 - 35.0/48.0) + 1.3125) - 35.0/96.0)
            - 0.7109375) - 7.0/768.0)*s);
    if (x < 2.5)
        return((x*(x*(x*(x*(x*(-1.0/20.0) + 7.0/12.0) - 2.625) + 133.0/24.0)
            - 5.140625) + 1267.0/960.0)*s);
    if (x < 3.5)
    {
        x *= 2.0;
        x -= 7.0;
        y = x*x;
        return((1.0/3840.0)*y*y*x*s);
    }
    return(0.0);
}

static double dwt7(double x)
{
    double y;
    int s;
    s = (x>0 ? 1 : -1);
    x = fabs(x);

    if (x < 1.0)
    {
        y = x*x;
        return(x*(y*(y*(x*(7.0/144.0) - 1.0/6.0) + 4.0/9.0) - 2.0/3.0)*s);
    }
    if (x < 2.0)
        return((x*(x*(x*(x*(x*(x*(-7.0/240.0) + 3.0/10.0)
            - 7.0/6.0) + 2.0) - 7.0/6.0) - 1.0/5.0) - 7.0/90.0)*s);
    if (x < 3.0)
        return((x*(x*(x*(x*(x*(x*(7.0/720.0) - 1.0/6.0)
            + 7.0/6.0) -38.0/9.0) + 49.0/6.0) - 23.0/3.0) + 217.0/90.0)*s);
    if (x < 4.0)
    {
        x -= 4;
        x *= x*x;
        x *= x;
        return((-1.0/720.0)*x*s);
    }
    return(0.0);
}

/***************************************************************************************
Generate B-spline basis functions
    d   - degree of spline
    x   - position relative to centre
    i   - pointer to first voxel position in convolution
    w   - vector of spline values

    Should really combine this function with wt2 to wt7 for most
    efficiency (as for case 0).

    Note that 0th degree B-spline returns nearest neighbour basis.
*/
static void weights(int d, double x, int *i, double w[])
{
    int k;

    *i = floor(x-(d-1)*0.5);
    x -= *i;

    switch (d){
    case 2:
        for(k=0; k<=2; k++) w[k] = wt2(x-k);
        break;
    case 3:
        for(k=0; k<=3; k++) w[k] = wt3(x-k);
        break;
    case 4:
        for(k=0; k<=4; k++) w[k] = wt4(x-k);
        break;
    case 5:
        for(k=0; k<=5; k++) w[k] = wt5(x-k);
        break;
    case 6:
        for(k=0; k<=6; k++) w[k] = wt6(x-k);
        break;
    case 7:
        for(k=0; k<=7; k++) w[k] = wt7(x-k);
        break;

    case 1:
        w[0] = 1.0-x;
        w[1] = x;
        break;
    case 0:
        w[0] = 1.0; /* Not correct at discontinuities */
        break;

    default:
        for(k=0; k<=7; k++) w[k] = wt7(x-k);
    }
}


/***************************************************************************************
Generate derivatives of B-spline basis functions
    d   - degree of spline
    x   - position relative to centre
    i   - pointer to first voxel position in convolution
    w   - vector of spline values

    Should really combine this function with dwt2 to dwt7 for most
    efficiency (as for case 0 and case 1).

    Note that 0th and 1st degree B-spline return derivatives of
    nearest neighbour and linear interpolation bases.
*/
static void dweights(int d, double x, int *i, double w[])
{
    int k;
    *i = floor(x-(d-1)*0.5);
    x -= *i;

    switch (d){
    case 2:
        for(k=0; k<=2; k++) w[k] = dwt2(x-k);
        break;
    case 3:
        for(k=0; k<=3; k++) w[k] = dwt3(x-k);
        break;
    case 4:
        for(k=0; k<=4; k++) w[k] = dwt4(x-k);
        break;
    case 5:
        for(k=0; k<=5; k++) w[k] = dwt5(x-k);
        break;
    case 6:
        for(k=0; k<=6; k++) w[k] = dwt6(x-k);
        break;
    case 7:
        for(k=0; k<=7; k++) w[k] = dwt7(x-k);
        break;

    case 1:
        w[0] = -1.0; /* Not correct at discontinuities */
        w[1] =  1.0; /* Not correct at discontinuities */
        break;
    case 0:
        w[0] = 0.0; /* Not correct at discontinuities */
        break;

    default:
        for(k=0; k<=7; k++) w[k] = dwt7(x-k);
    }
}


/***************************************************************************************
Work out what to do with positions outside the FOV
    i   - Co-ordinate (0<=i<m)
    m   - dimension
    returns reflected co-ordinate
*/
int mirror(int i, int m)
{
    int m2;
    i  = abs(i);
    if (i< m) return(i);
    if (m==1) return(0);
    m2 = (m-1)*2;
    i %= m2;
    return((i<m) ? i : m2-i);
}

/***************************************************************************************
Work out what to do with positions outside the FOV
        i       - Co-ordinate (0<=i<m)
        m       - dimension
        returns wrapped co-ordinate

        For MRI, it may be better to wrap the boundaries
        - especially in the read and phase encode directions.
*/
int wrap(int i, int m)
{
    if (i<0) return(m-1-((-i-1) % m));
    return(i % m);
}

/***************************************************************************************
Resample a point
    c   - Volume of B-spline coefficients
    m0,m1,m2    - dimensions of c
    x0,x1,x2    - co-ordinate to sample
    d   - degrees of splines used
    returns value of sampled point
*/
IMAGE_DTYPE sample3(IMAGE_DTYPE c[], int m0, int m1, int m2,
    IMAGE_DTYPE x0, IMAGE_DTYPE x1, IMAGE_DTYPE x2, int d[],
    int (*bnd[])())
{
    double w0[32], w1[32], w2[32]; /* B-spline weights */
    int    o0[32], o1[32], o2[32]; /* Offsets */
    int    i0,     i1,     i2;     /* Initial offsets */
    double d0,     d1,     d2;     /* Used by seperable convolution */
    int k;
    IMAGE_DTYPE *cp;

    /* Generate seperable B-spline basis functions */
    weights(d[0], x0, &i0, w0);
    weights(d[1], x1, &i1, w1);
    weights(d[2], x2, &i2, w2);

    /* Create lookups of voxel locations - for coping with edges */
    for(k=0; k<=d[0]; k++) o0[k] = bnd[0](k+i0, m0);
    for(k=0; k<=d[1]; k++) o1[k] = bnd[1](k+i1, m1)*m0;
    for(k=0; k<=d[2]; k++) o2[k] = bnd[2](k+i2, m2)*(m0*m1);

    /* Convolve coefficients with basis functions */
    d2 = 0.0;
    for(i2=0; i2<=d[2]; i2++)
    {
        d1 = 0.0;
        for(i1=0; i1<=d[1]; i1++)
        {
            cp = c+o2[i2]+o1[i1];
            d0 = 0.0;
            for(i0=0; i0<=d[0]; i0++)
                d0 += cp[o0[i0]] * w0[i0];
            d1 += d0 * w1[i1];
        }
        d2 += d1 * w2[i2];
    }
    return((IMAGE_DTYPE)d2);
}


/***************************************************************************************
Resample a point and its gradients
    c   - Volume of B-spline coefficients
    m0,m1,m2    - dimensions of c
    x0,x1,x2    - co-ordinate to sample
    d   - degrees of splines used
    pg0,pg1,pg2 - gradients
    returns value of sampled point
*/
IMAGE_DTYPE dsample3(IMAGE_DTYPE c[], int m0, int m1, int m2,
    IMAGE_DTYPE x0, IMAGE_DTYPE x1, IMAGE_DTYPE x2,
    int d[], IMAGE_DTYPE *pg0, IMAGE_DTYPE *pg1, IMAGE_DTYPE *pg2,
    int (*bnd[])())
{
    double  w0[32],  w1[32],  w2[32]; /* B-spline weights */
    double dw0[32], dw1[32], dw2[32]; /* B-spline derivatives */
    int     o0[32],  o1[32],  o2[32]; /* Offsets */
    int     i0,      i1,      i2;     /* Initial offsets */
    double  d0,      d1,      d2;     /* Used by seperable convolution */
    double g00, g10,g11, g20,g21,g22; /* Used for generating gradients */
    int k;
    IMAGE_DTYPE *cp;

    /* Generate seperable B-spline basis functions */
    weights(d[0], x0, &i0, w0);
    weights(d[1], x1, &i1, w1);
    weights(d[2], x2, &i2, w2);

    dweights(d[0], x0, &i0, dw0);
    dweights(d[1], x1, &i1, dw1);
    dweights(d[2], x2, &i2, dw2);

    /* Create lookups of voxel locations - for coping with edges */
    for(k=0; k<=d[0]; k++) o0[k] = bnd[0](k+i0, m0);
    for(k=0; k<=d[1]; k++) o1[k] = bnd[1](k+i1, m1)*m0;
    for(k=0; k<=d[2]; k++) o2[k] = bnd[2](k+i2, m2)*(m0*m1);

    /* Convolve coefficients with basis functions */
    g20 = g21 = g22 = d2 = 0.0;
    for(i2=0; i2<=d[2]; i2++)
    {
        g10 = g11 = d1 = 0.0;
        for(i1=0; i1<=d[1]; i1++)
        {
            cp = c+o2[i2]+o1[i1];
            g00 = d0  = 0.0;
            for(i0=0; i0<=d[0]; i0++)
            {
                d0  += cp[o0[i0]] *  w0[i0];
                g00 += cp[o0[i0]] * dw0[i0];
            }
            d1  += d0  *  w1[i1];
            g10 += g00 *  w1[i1];
            g11 += d0  * dw1[i1];
        }
        d2  += d1  *  w2[i2];
        g20 += g10 *  w2[i2];
        g21 += g11 *  w2[i2];
        g22 += d1  * dw2[i2];
    }
    *pg0 = g20;
    *pg1 = g21;
    *pg2 = g22;

    return((IMAGE_DTYPE)d2);
}

