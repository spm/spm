/*
 * $Id: spm_vol_utils.c 4452 2011-09-02 10:45:26Z guillaume $
 * John Ashburner
 */

#define TINY 5e-2

#ifdef SPM_UNSIGNED_CHAR
#define RESAMPLE resample_uchar
#define RESAMPLE_D resample_d_uchar
#define SLICE slice_uchar
#define RESAMPLE_0 resample_uchar_0
#define RESAMPLE_1 resample_uchar_1
#define RESAMPLE_D_1 resample_d_uchar_1
#define RESAMPLE_POLY resample_uchar_poly
#define RESAMPLE_D_POLY resample_d_uchar_poly
#define SLICE_0 slice_uchar_0
#define SLICE_1 slice_uchar_1
#define SLICE_POLY slice_uchar_poly
#define PLANE plane_uchar
#define GET(x) (x)
#define IMAGE_DTYPE unsigned char
#endif

#ifdef SPM_SIGNED_CHAR
#define RESAMPLE resample_schar
#define RESAMPLE_D resample_d_schar
#define SLICE slice_schar
#define RESAMPLE_0 resample_schar_0
#define RESAMPLE_1 resample_schar_1
#define RESAMPLE_D_1 resample_d_schar_1
#define RESAMPLE_POLY resample_schar_poly
#define RESAMPLE_D_POLY resample_d_schar_poly
#define SLICE_0 slice_schar_0
#define SLICE_1 slice_schar_1
#define SLICE_POLY slice_schar_poly
#define PLANE plane_schar
#define GET(x) (x)
#define IMAGE_DTYPE signed char
#endif

#ifdef SPM_SIGNED_SHORT
#ifdef SPM_BYTESWAP
#define GET(x) getshort(x)
#define RESAMPLE resample_short_s
#define RESAMPLE_D resample_d_short_s
#define SLICE slice_short_s
#define RESAMPLE_0 resample_short_s_0
#define RESAMPLE_1 resample_short_s_1
#define RESAMPLE_D_1 resample_d_short_s_1
#define RESAMPLE_POLY resample_short_s_poly
#define RESAMPLE_D_POLY resample_d_short_s_poly
#define SLICE_0 slice_short_s_0
#define SLICE_1 slice_short_s_1
#define SLICE_POLY slice_short_s_poly
#define PLANE plane_short_s
#else
#define GET(x) (x)
#define RESAMPLE resample_short
#define RESAMPLE_D resample_d_short
#define SLICE slice_short
#define RESAMPLE_0 resample_short_0
#define RESAMPLE_1 resample_short_1
#define RESAMPLE_D_1 resample_d_short_1
#define RESAMPLE_POLY resample_short_poly
#define RESAMPLE_D_POLY resample_d_short_poly
#define SLICE_0 slice_short_0
#define SLICE_1 slice_short_1
#define SLICE_POLY slice_short_poly
#define PLANE plane_short
#endif
#define IMAGE_DTYPE short int
#endif

#ifdef SPM_UNSIGNED_SHORT
#ifdef SPM_BYTESWAP
#define GET(x) getushort(x)
#define RESAMPLE resample_ushort_s
#define RESAMPLE_D resample_d_ushort_s
#define SLICE slice_ushort_s
#define RESAMPLE_0 resample_ushort_s_0
#define RESAMPLE_1 resample_ushort_s_1
#define RESAMPLE_D_1 resample_d_ushort_s_1
#define RESAMPLE_POLY resample_ushort_s_poly
#define RESAMPLE_D_POLY resample_d_ushort_s_poly
#define SLICE_0 slice_ushort_s_0
#define SLICE_1 slice_ushort_s_1
#define SLICE_POLY slice_ushort_s_poly
#define PLANE plane_ushort_s
#else
#define GET(x) (x)
#define RESAMPLE resample_ushort
#define RESAMPLE_D resample_d_ushort
#define SLICE slice_ushort
#define RESAMPLE_0 resample_ushort_0
#define RESAMPLE_1 resample_ushort_1
#define RESAMPLE_D_1 resample_d_ushort_1
#define RESAMPLE_POLY resample_ushort_poly
#define RESAMPLE_D_POLY resample_d_ushort_poly
#define SLICE_0 slice_ushort_0
#define SLICE_1 slice_ushort_1
#define SLICE_POLY slice_ushort_poly
#define PLANE plane_ushort
#endif
#define IMAGE_DTYPE unsigned short int
#endif

#ifdef SPM_SIGNED_INT
#ifdef SPM_BYTESWAP
#define GET(x) getint(x)
#define RESAMPLE resample_int_s
#define RESAMPLE_D resample_d_int_s
#define SLICE slice_int_s
#define RESAMPLE_0 resample_int_s_0
#define RESAMPLE_1 resample_int_s_1
#define RESAMPLE_D_1 resample_d_int_s_1
#define RESAMPLE_POLY resample_int_s_poly
#define RESAMPLE_D_POLY resample_d_int_s_poly
#define SLICE_0 slice_int_s_0
#define SLICE_1 slice_int_s_1
#define SLICE_POLY slice_int_s_poly
#define PLANE plane_int_s
#else
#define GET(x) (x)
#define RESAMPLE resample_int
#define RESAMPLE_D resample_d_int
#define SLICE slice_int
#define RESAMPLE_0 resample_int_0
#define RESAMPLE_1 resample_int_1
#define RESAMPLE_D_1 resample_d_int_1
#define RESAMPLE_POLY resample_int_poly
#define RESAMPLE_D_POLY resample_d_int_poly
#define SLICE_0 slice_int_0
#define SLICE_1 slice_int_1
#define SLICE_POLY slice_int_poly
#define PLANE plane_int
#endif
#define IMAGE_DTYPE int
#endif

#ifdef SPM_UNSIGNED_INT
#ifdef SPM_BYTESWAP
#define GET(x) getuint(x)
#define RESAMPLE resample_uint_s
#define RESAMPLE_D resample_d_uint_s
#define SLICE slice_uint_s
#define RESAMPLE_0 resample_uint_s_0
#define RESAMPLE_1 resample_uint_s_1
#define RESAMPLE_D_1 resample_d_uint_s_1
#define RESAMPLE_POLY resample_uint_s_poly
#define RESAMPLE_D_POLY resample_d_uint_s_poly
#define SLICE_0 slice_uint_s_0
#define SLICE_1 slice_uint_s_1
#define SLICE_POLY slice_uint_s_poly
#define PLANE plane_uint_s
#else
#define GET(x) (x)
#define RESAMPLE resample_uint
#define RESAMPLE_D resample_d_uint
#define SLICE slice_uint
#define RESAMPLE_0 resample_uint_0
#define RESAMPLE_1 resample_uint_1
#define RESAMPLE_D_1 resample_d_uint_1
#define RESAMPLE_POLY resample_uint_poly
#define RESAMPLE_D_POLY resample_d_uint_poly
#define SLICE_0 slice_uint_0
#define SLICE_1 slice_uint_1
#define SLICE_POLY slice_uint_poly
#define PLANE plane_uint
#endif
#define IMAGE_DTYPE unsigned int
#endif

#ifdef SPM_FLOAT
#ifdef SPM_BYTESWAP
#define GET(x) getfloat(x)
#define RESAMPLE resample_float_s
#define RESAMPLE_D resample_d_float_s
#define SLICE slice_float_s
#define RESAMPLE_0 resample_float_s_0
#define RESAMPLE_1 resample_float_s_1
#define RESAMPLE_D_1 resample_d_float_s_1
#define RESAMPLE_POLY resample_float_s_poly
#define RESAMPLE_D_POLY resample_d_float_s_poly
#define SLICE_0 slice_float_s_0
#define SLICE_1 slice_float_s_1
#define SLICE_POLY slice_float_s_poly
#define PLANE plane_float_s
#else
#define GET(x) (x)
#define RESAMPLE resample_float
#define RESAMPLE_D resample_d_float
#define SLICE slice_float
#define RESAMPLE_0 resample_float_0
#define RESAMPLE_1 resample_float_1
#define RESAMPLE_D_1 resample_d_float_1
#define RESAMPLE_POLY resample_float_poly
#define RESAMPLE_D_POLY resample_d_float_poly
#define SLICE_0 slice_float_0
#define SLICE_1 slice_float_1
#define SLICE_POLY slice_float_poly
#define PLANE plane_float
#endif
#define IMAGE_DTYPE float
#endif

#ifdef SPM_DOUBLE
#ifdef SPM_BYTESWAP
#define GET(x) getdouble(x)
#define RESAMPLE resample_double_s
#define RESAMPLE_D resample_d_double_s
#define SLICE slice_double_s
#define RESAMPLE_0 resample_double_s_0
#define RESAMPLE_1 resample_double_s_1
#define RESAMPLE_D_1 resample_d_double_s_1
#define RESAMPLE_POLY resample_double_s_poly
#define RESAMPLE_D_POLY resample_d_double_s_poly
#define SLICE_0 slice_double_s_0
#define SLICE_1 slice_double_s_1
#define SLICE_POLY slice_double_s_poly
#define PLANE plane_double_s
#else
#define GET(x) (x)
#define RESAMPLE resample_double
#define RESAMPLE_D resample_d_double
#define SLICE slice_double
#define RESAMPLE_0 resample_double_0
#define RESAMPLE_1 resample_double_1
#define RESAMPLE_D_1 resample_d_double_1
#define RESAMPLE_POLY resample_double_poly
#define RESAMPLE_D_POLY resample_d_double_poly
#define SLICE_0 slice_double_0
#define SLICE_1 slice_double_1
#define SLICE_POLY slice_double_poly
#define PLANE plane_double
#endif
#define IMAGE_DTYPE double
#endif

#include <math.h>
#include <stdlib.h>
#define RINT(A) floor((A)+0.5)
#include "spm_make_lookup.h"
#include "spm_getdata.h"

static void (*make_lookup)() = make_lookup_poly, (*make_lookup_grad)() = make_lookup_poly_grad;

/* Zero order hold resampling - nearest neighbour */
void RESAMPLE_0(m,vol,out,x,y,z,xdim,ydim,zdim,background, scale,offset)
int m, xdim,ydim,zdim;
double out[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int i;
    for (i=0; i<m; i++)
    {
        int xcoord, ycoord, zcoord;
        xcoord = floor(x[i]-0.5);
        ycoord = floor(y[i]-0.5);
        zcoord = floor(z[i]-0.5);
        if (xcoord>=0 && xcoord<xdim && ycoord>=0 &&
            ycoord<ydim && zcoord>=0 && zcoord<zdim)
            out[i] = scale[zcoord]*GET(vol[zcoord][xcoord  + xdim*ycoord])+offset[zcoord];
        else out[i] = background;
    }
}


/* First order hold resampling - trilinear interpolation */
void RESAMPLE_1(m,vol,out,x,y,z,xdim,ydim,zdim,background, scale,offset)
int m, xdim,ydim,zdim;
double out[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int i;
    for (i=0; i<m; i++)
    {
        double xi,yi,zi;
        xi=x[i]-1.0;
        yi=y[i]-1.0;
        zi=z[i]-1.0;

        if (    zi>=-TINY && zi<zdim+TINY-1 &&
            yi>=-TINY && yi<ydim+TINY-1 &&
            xi>=-TINY && xi<xdim+TINY-1)
        {
            double k111,k112,k121,k122,k211,k212,k221,k222;
            double dx1, dx2, dy1, dy2, dz1, dz2;
            int off1, off2, offx, offy, offz, xcoord, ycoord, zcoord;

            xcoord = (int)floor(xi); dx1=xi-xcoord; dx2=1.0-dx1;
            ycoord = (int)floor(yi); dy1=yi-ycoord; dy2=1.0-dy1;
            zcoord = (int)floor(zi); dz1=zi-zcoord; dz2=1.0-dz1;

            xcoord = (xcoord < 0) ? ((offx=0),0) : ((xcoord>=xdim-1) ? ((offx=0),xdim-1) : ((offx=1   ),xcoord));
            ycoord = (ycoord < 0) ? ((offy=0),0) : ((ycoord>=ydim-1) ? ((offy=0),ydim-1) : ((offy=xdim),ycoord));
            zcoord = (zcoord < 0) ? ((offz=0),0) : ((zcoord>=zdim-1) ? ((offz=0),zdim-1) : ((offz=1   ),zcoord));

            off1 = xcoord  + xdim*ycoord;
            off2 = off1+offy;
            k222 = GET(vol[zcoord     ][off1]); k122 = GET(vol[zcoord     ][off1+offx]);
            k212 = GET(vol[zcoord     ][off2]); k112 = GET(vol[zcoord     ][off2+offx]);
            k221 = GET(vol[zcoord+offz][off1]); k121 = GET(vol[zcoord+offz][off1+offx]);
            k211 = GET(vol[zcoord+offz][off2]); k111 = GET(vol[zcoord+offz][off2+offx]);

            /* resampled pixel value (trilinear interpolation) */
            out[i] =  (((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1)*scale[zcoord     ] + offset[zcoord     ])*dz2
                + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1)*scale[zcoord+offz] + offset[zcoord+offz])*dz1;
        }
        else out[i] = background;

    }
}

/* First order hold resampling - trilinear interpolation */
void RESAMPLE_D_1(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim,background, scale,offset)
int m, xdim,ydim,zdim;
double out[],gradx[],grady[],gradz[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int i;
    for (i=0; i<m; i++)
    {
        double xi,yi,zi;
        xi=x[i]-1.0;
        yi=y[i]-1.0;
        zi=z[i]-1.0;
        if (zi>=-TINY && zi<zdim+TINY-1 &&
            yi>=-TINY && yi<ydim+TINY-1 &&
            xi>=-TINY && xi<xdim+TINY-1)
        {
            double k111,k112,k121,k122,k211,k212,k221,k222;
            double dx1, dx2, dy1, dy2, dz1, dz2;
            int off1, off2, offx, offy, offz, xcoord, ycoord, zcoord;

            xcoord = (int)floor(xi); dx1=xi-xcoord; dx2=1.0-dx1;
            ycoord = (int)floor(yi); dy1=yi-ycoord; dy2=1.0-dy1;
            zcoord = (int)floor(zi); dz1=zi-zcoord; dz2=1.0-dz1;

            xcoord = (xcoord < 0) ? ((offx=0),0) : ((xcoord>=xdim-1) ? ((offx=0),xdim-1) : ((offx=1   ),xcoord));
            ycoord = (ycoord < 0) ? ((offy=0),0) : ((ycoord>=ydim-1) ? ((offy=0),ydim-1) : ((offy=xdim),ycoord));
            zcoord = (zcoord < 0) ? ((offz=0),0) : ((zcoord>=zdim-1) ? ((offz=0),zdim-1) : ((offz=1   ),zcoord));

            off1 = xcoord  + xdim*ycoord;
            off2 = off1+offy;
            k222 = GET(vol[zcoord     ][off1]); k122 = GET(vol[zcoord     ][off1+offx]);
            k212 = GET(vol[zcoord     ][off2]); k112 = GET(vol[zcoord     ][off2+offx]);
            k221 = GET(vol[zcoord+offz][off1]); k121 = GET(vol[zcoord+offz][off1+offx]);
            k211 = GET(vol[zcoord+offz][off2]); k111 = GET(vol[zcoord+offz][off2+offx]);

            /* resampled pixel value (trilinear interpolation) and gradients
               old code:
            gradx[i] = scale*(((k111     - k211    )*dy1 + (k121     - k221    )*dy2)*dz1
                            + ((k112     - k212    )*dy1 + (k122     - k222    )*dy2)*dz2);
            grady[i] = scale*(((k111*dx1 + k211*dx2)     - (k121*dx1 + k221*dx2)    )*dz1
                            + ((k112*dx1 + k212*dx2)     - (k122*dx1 + k222*dx2)    )*dz2);
            gradz[i] = scale*(((k111*dx1 + k211*dx2)*dy1 + (k121*dx1 + k221*dx2)*dy2)
                            - ((k112*dx1 + k212*dx2)*dy1 + (k122*dx1 + k222*dx2)*dy2)    );

            out[i]   = scale*(((k111*dx1 + k211*dx2)*dy1 + (k121*dx1 + k221*dx2)*dy2)*dz1
                            + ((k112*dx1 + k212*dx2)*dy1 + (k122*dx1 + k222*dx2)*dy2)*dz2) + offset; */

            gradx[i] = (((k111 - k211)*dy1 + (k121 - k221)*dy2)*scale[zcoord+offz])*dz1
                     + (((k112 - k212)*dy1 + (k122 - k222)*dy2)*scale[zcoord     ])*dz2;

            k111 = (k111*dx1 + k211*dx2)*scale[zcoord+offz]+offset[zcoord+offz];
            k121 = (k121*dx1 + k221*dx2)*scale[zcoord+offz]+offset[zcoord+offz];
            k112 = (k112*dx1 + k212*dx2)*scale[zcoord     ]+offset[zcoord     ];
            k122 = (k122*dx1 + k222*dx2)*scale[zcoord     ]+offset[zcoord     ];

            grady[i] = (k111 - k121)*dz1 + (k112 - k122)*dz2;

            k111 = k111*dy1 + k121*dy2;
            k112 = k112*dy1 + k122*dy2;

            gradz[i] = k111 - k112;
            out[i]   = k111*dz1 + k112*dz2;
        }
        else
        {
            out[i]   = background;
            gradx[i] = 0.0;
            grady[i] = 0.0;
            gradz[i] = 0.0;
        }
    }
}



/* Sinc resampling */
void RESAMPLE_POLY(m,vol,out,x,y,z,xdim,ydim,zdim, q,background, scale,offset)
int m, xdim,ydim,zdim, q;
double out[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int i;
    int dx1, dy1, dz1;
    static double tablex[255], tabley[255], tablez[255];

    for (i=0; i<m; i++)
    {
        if (z[i]>=1-TINY && z[i]<=zdim+TINY &&
            y[i]>=1-TINY && y[i]<=ydim+TINY &&
            x[i]>=1-TINY && x[i]<=xdim+TINY)
        {
            double dat=0.0, *tp1, *tp1end, *tp2end, *tp3end;

            make_lookup(x[i], q, xdim, &dx1, tablex, &tp3end);
            make_lookup(y[i], q, ydim, &dy1, tabley, &tp2end);
            make_lookup(z[i], q, zdim, &dz1, tablez, &tp1end);

            tp1 = tablez;
            dy1 *= xdim;

            while(tp1 <= tp1end)
            {
                IMAGE_DTYPE *dp2 = &vol[dz1][dy1];
                double dat2 = 0.0,
                *tp2 = tabley;
                while (tp2 <= tp2end)
                {
                    register double dat3 = 0.0, *tp3 = tablex;
                    register IMAGE_DTYPE *dp3 = dp2 + dx1;
                    while(tp3 <= tp3end)
                        dat3 += GET(*(dp3++)) * *(tp3++);
                    dat2 += dat3 * *(tp2++);
                    dp2  += xdim;
                }
                dat += (dat2*scale[dz1]+offset[dz1]) * *(tp1++);
                dz1 ++;
            }
            out[i] = dat;
        }
        else out[i] = background;
    }
}

/* Sinc resampling */
void RESAMPLE_D_POLY(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, q,background, scale,offset)
int m, xdim,ydim,zdim, q;
double out[],gradx[],grady[],gradz[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int i;
    int dx1, dy1, dz1;
    static double  tablex[255],  tabley[255],  tablez[255];
    static double dtablex[255], dtabley[255], dtablez[255];

    for (i=0; i<m; i++)
    {
        if (z[i]>=1-TINY && z[i]<=zdim+TINY &&
            y[i]>=1-TINY && y[i]<=ydim+TINY &&
            x[i]>=1-TINY && x[i]<=xdim+TINY)
        {
            double dat=0.0, datx = 0.0, daty = 0.0, datz = 0.0,
                  *tp1, *tp1end, *tp2end, *tp3end, *dp1;
            make_lookup_grad(x[i], q, xdim, &dx1, tablex, dtablex, &tp3end);
            make_lookup_grad(y[i], q, ydim, &dy1, tabley, dtabley, &tp2end);
            make_lookup_grad(z[i], q, zdim, &dz1, tablez, dtablez, &tp1end);
            tp1 =  tablez;
            dp1 = dtablez;
            dy1 *= xdim;
            while(tp1 <= tp1end)
            {
                IMAGE_DTYPE *d2 = &vol[dz1][dy1];
                double dat2  = 0.0, *tp2 =  tabley;
                double dat2x = 0.0, *dp2 = dtabley, dat2y = 0.0;
                while (tp2 <= tp2end)
                {
                    register IMAGE_DTYPE *d3 = d2 + dx1;
                    register double dat3  = 0.0, *tp3 =  tablex;
                    register double dat3x = 0.0, *dp3 = dtablex;
                    while(tp3 <= tp3end)
                    {
                        dat3x += GET(*(d3  )) * *(dp3++);
                        dat3  += GET(*(d3++)) * *(tp3++);
                    }
                    dat2x += dat3x * *(tp2  );
                    dat2  += dat3  * *(tp2++);
                    dat2y += dat3  * *(dp2++);
                    d2 += xdim;
                }

                datx += (dat2x*scale[dz1]) * *(tp1  );
                daty += (dat2y*scale[dz1]) * *(tp1  );
                dat2  = dat2*scale[dz1]+offset[dz1];
                dat  += dat2  * *(tp1++);
                datz += dat2  * *(dp1++);
                dz1  ++;
            }
            out[i]   = dat;
            gradx[i] = datx;
            grady[i] = daty;
            gradz[i] = datz;
        }
        else
        {
            out[i]   = background;
            gradx[i] = 0.0;
            grady[i] = 0.0;
            gradz[i] = 0.0;
        }
    }
}


/* Zero order hold resampling - nearest neighbour */
int SLICE_0(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale,offset)
double  mat[16], background, scale[],offset[];
double image[];
IMAGE_DTYPE *vol[];
int xdim1, ydim1, xdim2, ydim2, zdim2;
{
    double y, x2, y2, z2, s2, dx3=mat[0], dy3=mat[1], dz3=mat[2], ds3=mat[3];
    int t = 0;

    x2 = mat[12] + 0*mat[8];
    y2 = mat[13] + 0*mat[9];
    z2 = mat[14] + 0*mat[10];
    s2 = mat[15] + 0*mat[11];

    for(y=1; y<=ydim1; y++)
    {
        double x;
        double x3 = x2 + y*mat[4];
        double y3 = y2 + y*mat[5];
        double z3 = z2 + y*mat[6];
        double s3 = s2 + y*mat[7];

        for(x=1; x<=xdim1; x++)
        {
            int ix4, iy4, iz4;
            s3 += ds3;
            if (s3 == 0.0) return(-1);
            ix4 = floor(((x3 += dx3)/s3)-0.5);
            iy4 = floor(((y3 += dy3)/s3)-0.5);
            iz4 = floor(((z3 += dz3)/s3)-0.5);
            if (iz4>=0 && iz4<zdim2 && iy4>=0 && iy4<ydim2 && ix4>=0 && ix4<xdim2)
            {
                image[t] = scale[iz4]*(double)GET(vol[iz4][ix4 + xdim2*iy4]) + offset[iz4];
            }
            else image[t] = background;
            t++;
        }
    }
    return(0);
}

#define TINY 5e-2

/* First order hold resampling - trilinear interpolation */
int SLICE_1(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale,offset)
double  mat[16], background, scale[],offset[];
double image[];
IMAGE_DTYPE *vol[];
int xdim1, ydim1, xdim2, ydim2, zdim2;
{
    double y, x2, y2, z2, s2, dx3=mat[0], dy3=mat[1], dz3=mat[2], ds3=mat[3];
    int t = 0;

    x2 = mat[12] + 0*mat[8];
    y2 = mat[13] + 0*mat[9];
    z2 = mat[14] + 0*mat[10];
    s2 = mat[15] + 0*mat[11];

    for(y=1; y<=ydim1; y++)
    {
        double x;
        double x3 = x2 + y*mat[4];
        double y3 = y2 + y*mat[5];
        double z3 = z2 + y*mat[6];
        double s3 = s2 + y*mat[7];
        for(x=1; x<=xdim1; x++)
        {
            double x4,y4,z4;
            s3 += ds3;
            if (s3 == 0.0) return(-1);
            x4=(x3 += dx3)/s3-1.0;
            y4=(y3 += dy3)/s3-1.0;
            z4=(z3 += dz3)/s3-1.0;

            if (    z4>=-TINY && z4<zdim2+TINY-1 &&
                y4>=-TINY && y4<ydim2+TINY-1 &&
                x4>=-TINY && x4<xdim2+TINY-1)
            {
                double k111,k112,k121,k122,k211,k212,k221,k222;
                double dx1, dx2, dy1, dy2, dz1, dz2;
                int off1, off2, offx, offy, offz, ix4, iy4, iz4;

                ix4 = floor(x4); dx1=x4-ix4; dx2=1.0-dx1;
                iy4 = floor(y4); dy1=y4-iy4; dy2=1.0-dy1;
                iz4 = floor(z4); dz1=z4-iz4; dz2=1.0-dz1;

                ix4 = (ix4 < 0) ? ((offx=0),0) : ((ix4>=xdim2-1) ? ((offx=0),xdim2-1) : ((offx=1    ),ix4));
                iy4 = (iy4 < 0) ? ((offy=0),0) : ((iy4>=ydim2-1) ? ((offy=0),ydim2-1) : ((offy=xdim2),iy4));
                iz4 = (iz4 < 0) ? ((offz=0),0) : ((iz4>=zdim2-1) ? ((offz=0),zdim2-1) : ((offz=1    ),iz4));

                off1 = ix4  + xdim2*iy4;
                off2 = off1+offy;
                k222 = GET(vol[iz4     ][off1]); k122 = GET(vol[iz4     ][off1+offx]);
                k212 = GET(vol[iz4     ][off2]); k112 = GET(vol[iz4     ][off2+offx]);
                k221 = GET(vol[iz4+offz][off1]); k121 = GET(vol[iz4+offz][off1+offx]);
                k211 = GET(vol[iz4+offz][off2]); k111 = GET(vol[iz4+offz][off2+offx]);

                /* resampled pixel value (trilinear interpolation) */
                image[t] = (((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1)*scale[iz4     ] + offset[iz4     ])*dz2
                     + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1)*scale[iz4+offz] + offset[iz4+offz])*dz1;

            }
            else image[t] = background; 
            t++;
        }
    }
    return(0);
}


/* Sinc resampling */
int SLICE_POLY(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, q, background, scale,offset)
int ydim1,xdim1, xdim2,ydim2,zdim2, q;
double image[], mat[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int dx1, dy1, dz1;
    static double tablex[255], tabley[255], tablez[255];
    double y, x2, y2, z2, s2, dx3=mat[0], dy3=mat[1], dz3=mat[2], ds3=mat[3];

    x2 = mat[12] + 0*mat[8];
    y2 = mat[13] + 0*mat[9];
    z2 = mat[14] + 0*mat[10];
    s2 = mat[15] + 0*mat[11];

    for(y=1; y<=ydim1; y++)
    {
        double x;
        double x3 = x2 + y*mat[4];
        double y3 = y2 + y*mat[5];
        double z3 = z2 + y*mat[6];
        double s3 = s2 + y*mat[7];

        for(x=1; x<=xdim1; x++)
        {
            double x4,y4,z4;
            s3 += ds3;
            if (s3 == 0.0) return(-1);
            x4=(x3 += dx3)/s3;
            y4=(y3 += dy3)/s3;
            z4=(z3 += dz3)/s3;

            if (z4>=1-TINY && z4<=zdim2+TINY &&
                y4>=1-TINY && y4<=ydim2+TINY &&
                x4>=1-TINY && x4<=xdim2+TINY)
            {
                double dat=0.0, *tp1, *tp1end, *tp2end, *tp3end;

                make_lookup(x4, q, xdim2, &dx1, tablex, &tp3end);
                make_lookup(y4, q, ydim2, &dy1, tabley, &tp2end);
                make_lookup(z4, q, zdim2, &dz1, tablez, &tp1end);

                tp1 = tablez;
                dy1 *= xdim2;

                while(tp1 <= tp1end)
                {
                    IMAGE_DTYPE *dp2 = &vol[dz1][dy1];
                    double dat2 = 0.0,
                    *tp2 = tabley;
                    while (tp2 <= tp2end)
                    {
                        register double dat3 = 0.0, *tp3 = tablex;
                        register IMAGE_DTYPE *dp3 = dp2 + dx1;
                        while(tp3 <= tp3end)
                            dat3 += GET(*(dp3++)) * *(tp3++);
                        dat2 += dat3 * *(tp2++);
                        dp2 += xdim2;
                    }
                    dat += (dat2*scale[dz1]+offset[dz1]) * *(tp1++);
                    dz1 ++;
                }
                *(image++) = dat;
            }
            else *(image++) = background;
        }
    }
    return(0);
}

/* simple extraction of transverse plane */
void PLANE(p,image,vol,xdim,ydim,scale,offset)
int p, xdim,ydim;
double image[], scale[],offset[];
IMAGE_DTYPE *vol[];
{
    int n = xdim*ydim, i;
    IMAGE_DTYPE *ptr = vol[p-1];
    for(i=0; i<n; i++)
        image[i] = GET(ptr[i])*scale[p-1] + offset[p-1];
}

/* Extract a slice through the image */
int SLICE(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset)
int ydim1,xdim1, xdim2,ydim2,zdim2, hold;
double image[], mat[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    double t = 1e-10;
    /* attempt a nice easy transverse plane first */
    if (mat[0+0*4] > 1-t && mat[0+0*4] < 1+t &&
        mat[0+1*4] > -t  && mat[0+1*4] < t   &&
        mat[0+2*4] > -t  && mat[0+2*4] < t   &&
        mat[0+3*4] > -t  && mat[0+3*4] < t   &&
        mat[1+0*4] > -t  && mat[1+0*4] < t   &&
        mat[1+1*4] > 1-t && mat[1+1*4] < 1+t &&
        mat[1+2*4] > -t  && mat[1+2*4] < t   &&
        mat[1+3*4] > -t  && mat[1+3*4] < t   &&
        mat[2+0*4] > -t  && mat[2+0*4] < t   &&
        mat[2+1*4] > -t  && mat[2+1*4] < t   &&
        mat[2+2*4] > 1-t && mat[2+2*4] < 1+t &&
        fabs(RINT(mat[2+3*4])-mat[2+3*4])<t  &&
        mat[3+0*4] > -t  && mat[3+0*4] < t   &&
        mat[3+1*4] > -t  && mat[3+1*4] < t   &&
        mat[3+2*4] > -t  && mat[3+2*4] < t   &&
        mat[3+3*4] > 1-t && mat[3+3*4] < 1+t &&
        xdim1 == xdim2 && ydim1 == ydim2 && mat[2+3*4]>=1.0 && mat[2+3*4]<=zdim2)
    {
        int p;
        p = RINT(mat[2+3*4]);
        PLANE(p,image,vol,xdim2,ydim2,scale,offset);
        return(0);
    }
    else
    {
        if (hold<0)
        {
            hold=abs(hold);
            make_lookup = make_lookup_sinc;
            return(SLICE_POLY(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, hold+1, background, scale,offset));
        }
        else
        {
            make_lookup = make_lookup_poly;
        }
        if (hold == 0)
            return(SLICE_0(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale,offset));
        if (hold == 1)
            return(SLICE_1(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale,offset));
        else
            return(SLICE_POLY(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, hold+1, background, scale,offset));
    }
}

/* Resample image */
void RESAMPLE(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset)
int m, xdim,ydim,zdim, hold;
double out[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    if (hold<0)
    {
        hold=abs(hold);
        make_lookup = make_lookup_sinc;
        RESAMPLE_POLY(m,vol,out,x,y,z,xdim,ydim,zdim, hold+1, background, scale,offset);
    }
    else
    {
        make_lookup = make_lookup_poly;
        if (hold == 0)
            RESAMPLE_0(m,vol,out,x,y,z,xdim,ydim,zdim, background, scale,offset);
        else if (hold == 1)
            RESAMPLE_1(m,vol,out,x,y,z,xdim,ydim,zdim, background, scale,offset);
        else
            RESAMPLE_POLY(m,vol,out,x,y,z,xdim,ydim,zdim, hold+1, background, scale,offset);
    }
}

/* Resample image and derivatives */
void RESAMPLE_D(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset)
int m, xdim,ydim,zdim, hold;
double out[],gradx[],grady[],gradz[], x[], y[], z[], background, scale[],offset[];
IMAGE_DTYPE *vol[];
{
    if (hold<0)
    {
        hold=abs(hold);
        make_lookup_grad = make_lookup_sinc_grad;
        RESAMPLE_D_POLY(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold+1, background, scale,offset);
    }
    else
    {
        make_lookup_grad = make_lookup_poly_grad;
        if (hold == 1)
            RESAMPLE_D_1(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, background, scale,offset);
        else
            RESAMPLE_D_POLY(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold+1, background, scale,offset);
    }
}

