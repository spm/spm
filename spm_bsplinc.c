#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif
/*
 * This code is a modified version of that of Philippe Thevenaz, which I took from:
 *	http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs.
 *
 * See:
 *	M. Unser, A. Aldroubi and M. Eden.
 *	"B-Spline Signal Processing: Part I-Theory,"
 *	IEEE Transactions on Signal Processing 41(2):821-832 (1993).
 *
 *	M. Unser, A. Aldroubi and M. Eden.
 *	"B-Spline Signal Processing: Part II-Efficient Design and Applications,"
 *	IEEE Transactions on Signal Processing 41(2):834-848 (1993).
 *
 *	M. Unser.
 *	"Splines: A Perfect Fit for Signal and Image Processing,"
 *	IEEE Signal Processing Magazine 16(6):22-38 (1999).
 *
*/

#include <math.h>
#include <mex.h>
#include "spm_sys_deps.h"
#include "spm_mapping.h"
#include "spm_datatypes.h"


/***************************************************************************************
Starting mirrored boundary condition based on Eq. 2.6 of Unser's 2nd 1993 paper.
	c - vector of unfiltered data
	m - length of c
	p - pole (root of polynomial)
	function returns value that c[0] should initially take
*/
static double cc_mirror(double c[], int m, double p)
{
	double s, pi, p2i, ip;
	int    i, m1;

	/* Initial causal coefficient assuming reflected boundaries */
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
End mirrored boundary condition
	c - first pass filtered data
	m - length of filtered data (must be > 1)
	p - pole
	function returns value for c[m-1] before 2nf filter pass
*/
static double icc_mirror(double c[],int m, double p)
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
	/* compute gain */
	int j;
	double lambda = 1.0;
	for (j = 0; j < np; j++)
		lambda = lambda*(1.0-p[j])*(1.0-1.0/p[j]);
	return(lambda);
}


/***************************************************************************************
One dimensional recursive filtering - assuming mirror boundaries
See Eq. 2.5 of Unsers 2nd 1993 paper.
	c - original vector on input, coefficients on output
	m - length of vector
	p - poles (polynomial roots)
	np - number of poles
*/
static void splinc(double c[], int m, double p[], int np)
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
	 o - order of B-spline
	 np - number of roots of magnitude less than one
	 p - roots.
*/
static int get_poles(int o, int *np, double p[])
{
	/* Return polynomial roots that are less than one. */
	switch (o) {
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
			p[0] = -0.36134122590022010879;
			p[1] = -0.013725429297339124604;
			break;
		case 5: /* roots([1 26 66 26 1]) */
			*np   = 2;
			p[0] = -0.43057534709997380418;
			p[1] = -0.043096288203264644656;
			break;
		case 6: /* roots([1 722 10543 23548 10543 722 1]) */
			*np   = 3;
			p[0] = -0.4882945893030444795;
			p[1] = -0.081679271076237569549;
			p[2] = -0.0014141518083258179488;
			break;
		case 7: /* roots([1 120 1191 2416 1191 120 1]) */
			*np   = 3;
			p[0] = -0.53528043079643827795;
			p[1] = -0.12255461519232672962;
			p[2] = -0.0091486948096082768705;
			break;
		default:
			return(1);
	}
	return(0);
}


/***************************************************************************************
Deconvolve the B-spline basis functions from the image volume
	vol - a handle for the volume to deconvolve
	c - the coefficients (arising from the deconvolution)
	o - the spline order
*/
static int vol_coeffs(MAPTYPE *vol, double c[], int o)
{
	double	p[4], *cp;
	int	np;
	int	i, j, k, n;
	double f[10240];

	if (vol->dim[0]>10240 || vol->dim[1]>10240 ||vol->dim[2]>10240)
		return(1);

	if (get_poles(o, &np, p))
		return(1);

	if (o<=1) /* Just do a straight copy */
	{
		cp = c;
		for(k=0; k<vol->dim[2]; k++)
		{
			double dk = k+1;
			for(j=0; j<vol->dim[1]; j++)
			{
				double dj = j+1;
				for(i=0;i<vol->dim[0];i++, cp++)
				{
					double di = i+1;
					resample(1,vol,cp,&di,&dj,&dk,0, 0.0);
				}
			}
		}
	}
	else /* Deconvolve along the fastest dimension (X) */
	{
		for(k=0; k<vol->dim[2]; k++)
		{
			double dk = k+1;
			for(j=0; j<vol->dim[1]; j++)
			{
				double dj = j+1;
				for(i=0;i<vol->dim[0];i++)
				{
					double di = i+1;
					resample(1,vol,&f[i],&di,&dj,&dk,0, 0.0);
				}
				splinc(f, vol->dim[0], p, np);
				cp = &c[vol->dim[0]*(j+vol->dim[1]*k)];
				for(i=0;i<vol->dim[0];i++, cp++)
					*cp = f[i];
			}
		}
	}

	if (o>1) /* Deconvolve in middle dimension (Y) */
	{
		n =vol->dim[0];
		for(k=0; k<vol->dim[2]; k++)
		{
			for(i=0;i<vol->dim[0];i++)
			{
				cp = &c[i+vol->dim[0]*vol->dim[1]*k];
				for(j=0; j<vol->dim[1]; j++, cp+=n)
					f[j] = *cp;
				splinc(f, vol->dim[1], p, np);
				cp = &c[i+vol->dim[0]*vol->dim[1]*k];
				for(j=0; j<vol->dim[1]; j++, cp+=n)
					*cp = f[j];
			}
		}
	}

	if (o>1) /* Deconvolve in the slowest dimension (Z) */
	{
		n = vol->dim[0]*vol->dim[1];
		for(j=0; j<vol->dim[1]; j++)
		{
			for(i=0;i<vol->dim[0];i++)
			{
				cp = &c[i+vol->dim[0]*j];
				for(k=0; k<vol->dim[2]; k++, cp+=n)
					f[k] = *cp;
				splinc(f, vol->dim[2], p, np);
				cp = &c[i+vol->dim[0]*j];
				for(k=0; k<vol->dim[2]; k++, cp+=n)
					*cp = f[k];
			}
		}
	}

	return(0);
}

/***************************************************************************************
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int k, o;
	MAPTYPE *vol, *get_maps();
	double *c;

	if (nrhs < 2 || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");
	if (mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || mxGetM(prhs[1])*mxGetN(prhs[1]) != 1)
		mexErrMsgTxt("Inappropriate usage.");

	o = rint(mxGetPr(prhs[1])[0]);
	if (o<0 || o>7)
		mexErrMsgTxt("Bad spline order.");

	vol=get_maps(prhs[0], &k);
	if (k!=1)
	{
		free_maps(vol, k);
		mexErrMsgTxt("Too many images.");
	}

	plhs[0] = mxCreateNumericArray(3,vol->dim, mxDOUBLE_CLASS, mxREAL);
	c = mxGetPr(plhs[0]);

	if (vol_coeffs(vol, c, o))
	{
		free_maps(vol, k);
		mexErrMsgTxt("Problem with deconvolution.");
	}
	free_maps(vol, k);
}
