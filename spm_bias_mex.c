#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include "spm_mapping.h"
#include "spm_sys_deps.h"
#include <math.h>

#define EPS 2.2204460492503130808e-16
#define MAXB 128 /* Maximum number of bases functions in each dimension */

#define IGNORE_ZEROS


/********************************************************************************/
/* Extract a slice from a mapped volume
 * vol - mapped volume
 * i   - slice number (first slice is 0)
 * dat - extracted slice
 *
 * required: vol, i
 * modified: dat
*/
static void get_slice(MAPTYPE *vol, int i, double dat[])
{
	static double mat[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
	mat[14] = i+1.0;
	slice(mat,dat,vol->dim[0],vol->dim[1], vol, 0,0);

}
/********************************************************************************/
/* Normalise the histogram so that it integrates to unity, taking the width of
 * the bins into account.
 * nh       - number of histogram bins
 * h        - histogram
 * binwidth - width of each histogram bin (1/scale)
 *
 * required: nh, h, binwidth
 * modified: h
 */
static void normhist(int nh, double h[], double binwidth, int *psh)
{
	double sh = 0.0;
	int i;
	for(i=0; i<nh; i++)
		sh += h[i];
	for(i=0; i<nh; i++)
		h[i]/= binwidth*sh;
	*psh = rint(sh);
}

/********************************************************************************/
/* Create a histogram from data modulated by a linear combination of seperable
 * basis functions
 * m0  - number of bases in x dimension
 * m1  - number of bases in y dimension
 * m2  - number of bases in z dimension
 * n0  - x dimension of volume
 * n1  - y dimension of volume
 * n2  - z dimension of volume
 * vol - mapped volume
 * B0  - seperable bases - x
 * B1  - seperable bases - y
 * B2  - seperable bases - z
 * f   - basis function coefficients
 * nh  - number of histogram bins
 * s   - scale from image intensity to histogram bin number (1/binwidth)
 * h   - histogram
 * dat - buffer for storing a plane of data
 *
 * required: m0, m1, m2, n0, n1, n2, vol, B0, B1, B2, f, nh, s
 * modified: h, dat
*/
static void hist(int m0, int m1, int m2, int n0, int n1, int n2,
	MAPTYPE *vol, double *B0, double *B1, double *B2, double *f,
	int nh, double s, double *h, double *dat, int *psh)
{
	int i0, i1, i2;
	int j0, j1, j2;
	double t2[MAXB*MAXB], t1[MAXB];

	for(j2=0;j2<n2; j2++)
	{
		get_slice(vol,j2,dat);
		for(i0=0; i0<m0*m1; i0++)
		{
			t2[i0] = 0.0;
			for(i2=0;i2<m2; i2++)
				t2[i0] += B2[j2+n2*i2]*f[i0+m0*m1*i2];
		}
		for(j1=0; j1<n1; j1++)
		{
			double *pdat = dat+j1*n0;
			for(i0=0; i0<m0; i0++)
			{
				t1[i0] = 0.0;
				for(i1=0;i1<m1; i1++)
					t1[i0] += B1[j1+n1*i1]*t2[i0+m0*i1];
			}
			for(j0=0; j0<n0; j0++)
			{
				double y, y0 = pdat[j0]*s;
#ifdef IGNORE_ZEROS
				if (y0>0.0)
#else
				if (y0>=0.0)
#endif
				{
					int iy;
					double dy, sc = 0.0;
					for(i0=0;i0<m0; i0++)
						sc += B0[j0+n0*i0]*t1[i0];
					y = y0*exp(sc);
					if (y>=0 && y<(nh-1))
					{
						iy = floor(y);
						dy = y-iy;
						h[iy] += 1-dy;
						h[iy+1] += dy;
					}
				}
			}
		}
	}
	normhist(nh,h,1/s,psh);
}
/********************************************************************************/

/********************************************************************************/
/* Create voxel values in weighting images and increment likelihoods using
 * lookup tables
 * dat - raw voxel intensity
 * sc  - estimated scaling
 * nh  - number of bins in lookup table
 * s   - scale from image intensity to bin number
 * wt0 - first weight
 * wt1 - second weight
 * ll  - incremented log-likelihood
 *
 * required: dat, sc, nh, s
 * modified: wt0, wt1, ll
 */
static void weights(double dat, double sc,
	int nh, double s, double h[],
	double *wt0, double *wt1, double *ll)
{
	double y, dy, f, df, t;
	int iy;

	y = dat*sc*s;
#ifdef IGNORE_ZEROS
	if (y>0.0 && y<nh-1.0)
#else
	if (y>=0.0 && y<nh-1.0)
#endif
	{
		iy   = (int)floor(y);
		dy   = y-iy;
		f    = ((1.0-dy)*h[iy] + dy*h[iy+1])+EPS;
		df   = (h[iy+1]-h[iy])*s;
		t    = dat*sc*df/f;
		*wt0 = -1.0 - t;
		*wt1 =  t*t - t;
		*ll -=  log(f*sc);
	}
	else
	{
		/*
		*wt0 = -1.0/sc;
		*wt1 =  1.0/(sc*sc);
		*/
		*wt0 = 0.0;
		*wt1 = 0.0;
	}
}
/********************************************************************************/

/********************************************************************************/
/* Beta = kron(B1,B0)'*img(:)
 * m0  - rows in B0
 * m1  - rows in B1
 * n0  - columns in B0
 * n1  - columns in B1
 * img - m0*m1 vector
 * B0  - basis functions - x
 * B1  - basis functions - y
 * Beta - resulting vector
 *
 * required: n0, n1, m0, m1, img, B0, B1
 * modified: Beta
*/
static void kronutil1(int m0, int m1, int n0, int n1,
	double img[], double B0[], double B1[], double Beta[])
{
	int j1,j2, i1,i2;
	double Beta1[MAXB];

	/* Zero Beta */
	for(j1=0; j1<n0*n1; j1++) Beta[j1]=0.0;

	for(i2=0; i2<m1; i2++)
	{
		/* generate small Beta */
		for(j1=0; j1<n0; j1++) Beta1[j1]=0.0;
		for(i1=0; i1<m0; i1++)
		{
			double wt = img[i1+m0*i2];
			for(j1=0; j1<n0; j1++)
				Beta1[j1] += B0[i1+j1*m0]*wt;
		}

		/* kronecker tensor product to increment Beta */
		for(j2=0; j2<n1; j2++)
		{
			double wt = B1[i2+j2*m1];
			double *ptrb = Beta+(n0*j2);
			for(j1=0; j1<n0; j1++)
				ptrb[j1] += wt*Beta1[j1];
		}
	}
}
/********************************************************************************/

/********************************************************************************/
/* Alpha = kron(B1,B0)'*diag(img(:))*kron(B1,B0)
 * m0  - rows in B0
 * m1  - rows in B1
 * n0  - columns in B0
 * n1  - columns in B1
 * img - m0*m1 vector
 * B0  - basis functions - x
 * B1  - basis functions - y
 * Alpha - resulting matrix
 *
 * required: n0, n1, m0, m1, img, B0, B1
 * modified: Alpha
*/
static void kronutil2(int m0, int m1, int n0, int n1,
	double img[], double B0[], double B1[], double Alpha[])
{
	int j11,j12, j21,j22, i1, i2;
	double Alpha1[MAXB*MAXB];

	/* Zero Alpha */
	for(j21=0; j21<n0*n1; j21++)
		for(j11=0; j11<=j21; j11++)
			Alpha[j11+j21*n0*n1]=0.0;

	for(i2=0; i2<m1; i2++)
	{
		/* zero small Alpha */
		for(j21=0; j21<n0; j21++)
			for(j11=0; j11<=j21; j11++)
				Alpha1[j11+j21*n0]=0.0;

		/* generate upper half of small Alpha */
		for(i1=0; i1<m0; i1++)
		{
			double wt = img[i1+m0*i2];
			for(j21=0; j21<n0; j21++)
			{
				double wt2 = wt*B0[i1+j21*m0];
				for(j11=0; j11<=j21; j11++)
					Alpha1[j11+j21*n0] += wt2*B0[i1+j11*m0];
			}
		}

		/* kronecker tensor product to increment upper half of large Alpha */
		for(j22=0; j22<n1; j22++)
		{
			double wt = B1[i2+j22*m1];
			for(j12=0; j12<=j22; j12++)
			{
				double *ptra = Alpha+(n0*j12+n0*n0*n1*j22);
				double wt2 = wt*B1[i2+j12*m1];

				for(j21=0; j21<n0; j21++)
					for(j11=0; j11<=j21; j11++)
						ptra[j11+j21*n0*n1] += wt2*Alpha1[j11+j21*n0];
			}
		}
	}
	/* fill in lower triangle of symmetric submatrices of Alpha */
	for(j22=0; j22<n1; j22++)
		for(j12=0; j12<j22; j12++)
		{
			double *ptra = Alpha+(n0*j12+n0*n0*n1*j22);
			for(j21=0; j21<n0; j21++)
				for(j11=0; j11<j21; j11++)
					ptra[j21+j11*n0*n1] = ptra[j11+j21*n0*n1];
		}

	/* fill in lower triangle of Alpha */
	for(j22=0; j22<n1*n0; j22++)
		for(j12=0; j12<j22; j12++)
			Alpha[j22+j12*n0*n1]=Alpha[j12+j22*n0*n1];
}
/********************************************************************************/

/********************************************************************************/
/* Increments a Alpha and Beta from one plane of data
 * Beta  += kron(B2,kron(B1,B0))'*diag(wt0(:))
 * Alpha += kron(B2,kron(B1,B0))'*diag(wt1(:))*kron(B2,kron(B1,B0))
 * m0  - rows in B0
 * m1  - rows in B1
 * m2  - rows in B2
 * n0  - columns in B0
 * n1  - columns in B1
 * n2  - columns in B2
 * i2  - row of interest in B2
 * wt0 - m0*m1 vector
 * wt1 - m0*m1 vector
 * B0  - basis functions - x
 * B1  - basis functions - y
 * B2  - basis functions - z
 * Alpha1 - buffer
 * Beta1  - buffer
 * Alpha  - resulting matrix
 * Beta   - resulting vector
 *
 * required: m0, m1, m2, n0, n1, n2, i2, wt0, wt1, B0, B1, B2
 * modified: Alpha1, Beta1, Alpha, Beta
*/
static void kronutil(int m0, int m1, int m2, int n0, int n1, int n2, int i2,
	double wt0[], double wt1[], double B0[], double B1[], double B2[],
	double Alpha1[], double Beta1[], double Alpha[], double Beta[])
{
	int j01,j2, n01;

	kronutil1(m0, m1, n0, n1, wt0, B0, B1, Beta1);
	kronutil2(m0, m1, n0, n1, wt1, B0, B1, Alpha1);

	n01 = n0*n1;
	for(j2=0; j2<n2; j2++)
	{
		int k01,k2;
		for(j01=0; j01<n01; j01++)
			Beta[j01 + j2*n01] += Beta1[j01]*B2[i2 + j2*m2];

		for(k2=0; k2<n2; k2++)
		{
			double tmp = B2[i2 + j2*m2]*B2[i2 + k2*m2];
			for(j01=0; j01<n01; j01++)
				for(k01=0; k01<n01; k01++)
					Alpha[k01 + k2*n01 + (j01 + j2*n01)*(n01*n2)] += Alpha1[k01 + j01*n01]*tmp;
		}
	}
}
/********************************************************************************/

/********************************************************************************/
static double iteration(int m0, int m1, int m2, int n0, int n1, int n2,
	double Alpha[], double Beta[], MAPTYPE *vol,
	double B0[], double B1[], double B2[], double f[], double mx, int nh, double h[], int *psh)
{
	int j2, j1, j0, i2, i1, i0;
	double *dat, *wt0, *wt1, *Alpha1, *Beta1;
	double t2[MAXB*MAXB], t1[MAXB], ll;

	dat    = (double *)mxCalloc(n1*n0, sizeof(double));
	wt0    = (double *)mxCalloc(n1*n0, sizeof(double));
	wt1    = (double *)mxCalloc(n1*n0, sizeof(double));
	Alpha1 = (double *)mxCalloc(m1*m0*m1*m0, sizeof(double));
	Beta1  = (double *)mxCalloc(m1*m0, sizeof(double));

	hist(m0,m1,m2, n0,n1,n2, vol, B0,B1,B2,f, nh,(nh-1.0)/mx,h, dat, psh);

	ll = 0.0;

	for(j2=0;j2<n2; j2++)
	{
		get_slice(vol,j2,dat);
		for(i0=0; i0<m0*m1; i0++)
		{
			t2[i0] = 0.0;
			for(i2=0;i2<m2; i2++)
				t2[i0] += B2[j2+n2*i2]*f[i0+m0*m1*i2];
		}
		for(j1=0; j1<n1; j1++)
		{
			for(i0=0; i0<m0; i0++)
			{
				t1[i0] = 0.0;
				for(i1=0;i1<m1; i1++)
					t1[i0] += B1[j1+n1*i1]*t2[i0+m0*i1];
			}
			for(j0=0; j0<n0; j0++)
			{
				double sc = 0.0;
				for(i0=0; i0<m0; i0++)
					sc += B0[j0+n0*i0]*t1[i0];
				weights(dat[j0+n0*j1], exp(sc), nh, (nh-1.0)/mx,h, &wt0[j0+n0*j1], &wt1[j0+n0*j1], &ll);
			}
		}
		kronutil(n0,n1,n2, m0,m1,m2, j2, wt0,wt1, B0,B1,B2, Alpha1,Beta1, Alpha,Beta);
	}

	mxFree((void *)dat);
	mxFree((void *)wt0);
	mxFree((void *)wt1);
	mxFree((void *)Alpha1);
	mxFree((void *)Beta1);

	return(ll);
}
/********************************************************************************
* Unused code */
double constr(int m0, int m1, int m2, int n0, int n1, int n2,
	double Beta[], MAPTYPE *vol,
	double B0[], double B1[], double B2[])
{
	int j2, i2, i01;
	double *dat, *Beta1, sm;
	dat    = (double *)mxCalloc(n1*n0, sizeof(double));
	Beta1  = (double *)mxCalloc(m1*m0, sizeof(double));

	sm = 0.0;
	for(j2=0;j2<n2; j2++)
	{
		get_slice(vol,j2,dat);
		kronutil1(n0,n1,m0,m1, dat, B0, B1, Beta1);
		for(i2=0; i2<m2; i2++)
			for(i01=0; i01<m0*m1; i01++)
				Beta[i01 + i2*m0*m1] += Beta1[i01]*B2[j2 + i2*n2];
		for(i01=0; i01<n0*n1; i01++)
			sm += dat[i01];
	}
	mxFree((void *)Beta1);
	mxFree((void *)dat);
	return(sm);
}

void mexFunction1(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *vol, *get_maps();
	int m0,m1,m2, n0,n1,n2;
	double *B0, *B1, *B2;
	for(n0=1; n0<nrhs; n0++)
	{
		if (!mxIsNumeric(prhs[n0]) || mxIsComplex(prhs[n0]) ||
			mxIsSparse(prhs[n0]) || !mxIsDouble(prhs[n0]))
		{
			free_maps(vol, 1);
			mexErrMsgTxt("Arguments must be numeric, real, full and double.");
		}
	}

	vol = get_maps(prhs[0], &n0);
	if (n0!=1)
	{
		free_maps(vol, n0);
		mexErrMsgTxt("Inappropriate usage.");
	}
	n0 = vol->dim[0];
	n1 = vol->dim[1];
	n2 = vol->dim[2];

	m0 = mxGetN(prhs[1]);
	m1 = mxGetN(prhs[2]);
	m2 = mxGetN(prhs[3]);

	if (m0>MAXB || m1>MAXB || m2>MAXB)
	{
		free_maps(vol, 1);
		mexErrMsgTxt("Too many basis functions.");
	}

	B0 = mxGetPr(prhs[1]);
	B1 = mxGetPr(prhs[2]);
	B2 = mxGetPr(prhs[3]);

	if (mxGetM(prhs[1])!=n0 || mxGetM(prhs[2])!=n1 || mxGetM(prhs[3])!=n2)
	{
		free_maps(vol, 1);
		mexErrMsgTxt("Inappropriate matrix dimensions.");
	}
	plhs[0] = mxCreateDoubleMatrix(m0*m1*m2,1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
	mxGetPr(plhs[1])[0] = constr(m0,m1,m2,n0,n1,n2,mxGetPr(plhs[0]),vol,B0,B1,B2);
	free_maps(vol, 1);
}
/********************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m0, m1, m2,  n0, n1, n2, i, nh, sh;
	MAPTYPE *vol, *get_maps();
	double *B0, *B1, *B2, mx, *f, *Alpha, *Beta, *h;

	if (nrhs == 4 && nlhs <=2)
	{
		mexFunction1(nlhs, plhs, nrhs, prhs);
		return;
	}
	if (nrhs != 6 || nlhs > 5) mexErrMsgTxt("Inappropriate usage.");

	vol = get_maps(prhs[0], &n0);
	if (n0!=1)
	{
		free_maps(vol, n0);
		mexErrMsgTxt("Inappropriate usage.");
	}

	for(n0=1; n0<nrhs; n0++)
	{
		if (!mxIsNumeric(prhs[n0]) || mxIsComplex(prhs[n0]) ||
			mxIsSparse(prhs[n0]) || !mxIsDouble(prhs[n0]))
		{
			free_maps(vol, 1);
			mexErrMsgTxt("Arguments must be numeric, real, full and double.");
		}
	}

	n0 = vol->dim[0];
	n1 = vol->dim[1];
	n2 = vol->dim[2];

	m0 = mxGetN(prhs[1]);
	m1 = mxGetN(prhs[2]);
	m2 = mxGetN(prhs[3]);

	if (mxGetM(prhs[1])!=n0 || mxGetM(prhs[2])!=n1 || mxGetM(prhs[3])!=n2
		|| mxGetM(prhs[4])*mxGetN(prhs[4]) != m0*m1*m2
		|| mxGetM(prhs[5])*mxGetN(prhs[5])!=2)
	{
		free_maps(vol, 1);
		mexErrMsgTxt("Inappropriate matrix dimensions.");
	}

	if (m0>MAXB || m1>MAXB || m2>MAXB)
	{
		free_maps(vol, 1);
		mexErrMsgTxt("Too many basis functions.");
	}

	B0 = mxGetPr(prhs[1]);
	B1 = mxGetPr(prhs[2]);
	B2 = mxGetPr(prhs[3]);
	f  = mxGetPr(prhs[4]);
	mx = mxGetPr(prhs[5])[0];
	nh = floor(mxGetPr(prhs[5])[1]);
	if (nh<1 || nh > 10240)
	{
		free_maps(vol, 1);
		mexErrMsgTxt("Inappropriate number of bins in histogram.");
	}

	plhs[0] = mxCreateDoubleMatrix(m0*m1*m2,m0*m1*m2, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(m0*m1*m2,1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(nh,1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);

	/* Initialised to zero */
	Alpha   = mxGetPr(plhs[0]);
	Beta    = mxGetPr(plhs[1]);

	h = mxGetPr(plhs[3]);

	mxGetPr(plhs[2])[0] = iteration(m0,m1,m2, n0,n1,n2, Alpha,Beta, vol, B0,B1,B2,f,mx, nh,h, &sh);
	mxGetPr(plhs[4])[0] = (double)sh;
	free_maps(vol, 1);
}
