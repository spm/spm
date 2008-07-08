/*
 * $Id: spm_bias_mex.c 1893 2008-07-08 15:05:40Z john $
 * John Ashburner
 */

#include <math.h>
#include "mex.h"
#include "spm_mapping.h"
#include "spm_datatypes.h"
#define RINT(A) floor((A)+0.5)
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

	if ((vol->dtype == SPM_UNSIGNED_CHAR)    || (vol->dtype == SPM_SIGNED_SHORT)   ||
	    (vol->dtype == SPM_SIGNED_INT)       || (vol->dtype == SPM_SIGNED_SHORT_S) ||
	    (vol->dtype == SPM_SIGNED_INT_S)     || (vol->dtype == SPM_SIGNED_CHAR)    ||
	    (vol->dtype == SPM_UNSIGNED_SHORT)   || (vol->dtype == SPM_UNSIGNED_INT)   ||
	    (vol->dtype == SPM_UNSIGNED_SHORT_S) || (vol->dtype == SPM_UNSIGNED_INT_S))
	{
		/* Only a relatively small finite number of values can be stored as
		   integers.  This can cause aliasing problems in the histograms.
		   The ideal solution is to assume an intensity PDF for each voxel
		   that is a top-hat, where the width is the scalefactor.
		   Adding to the histogram should really use the linear interpolation
		   kernel convolved with a top-hat of the same width as the scalefactor,
		   Instead, I've added some pseudo-random noise - much easier to implement
		   slightly faster, but not as good */

		static float ran[] = {
			0.141571,-0.193651,0.160932,-0.141983,0.438199,-0.012333,-0.40901,0.173834,
			0.0148803,-0.278421,0.225009,-0.431753,0.464124,-0.292343,-0.338882,0.138221,
			-0.499772,-0.164367,-0.2249,-0.455472,-0.406104,-0.0900026,0.316892,0.370517,
			-0.477445,0.227177,0.348009,0.2286,0.455099,0.156351,0.242305,-0.155034,
			0.384022,-0.152756,-0.440523,0.218415,0.458214,-0.343165,-0.0836469,-0.405965,
			-0.0500535,0.369152,-0.108383,-0.247216,-0.145618,0.242978,0.150832,0.439793,
			0.332799,-0.0300221,0.129866,-0.441812,0.0421874,-0.044274,0.363087,0.355197,
			-0.0277443,0.286924,0.155982,-0.49996,-0.368763,-0.00512524,-0.461667,-0.272564,
			-0.172117,0.399469,-0.18627,-0.248324,-0.0670109,0.342382,-0.315511,0.00817921,
			-0.0477603,-0.174416,-0.119924,0.38648,0.261261,0.383766,-0.0425937,0.299202,
			-0.365923,-0.434686,-0.124855,-0.126477,-0.0159776,0.469459,-0.157939,-0.247311,
			0.0848869,0.0237036,-0.336581,-0.0136019,-0.00393925,0.343194,0.306198,0.357786,
			0.109754,0.0657304,0.111899,-0.397023,-0.341684,-0.0863501,0.0604104,-0.231323,
			0.284254,-0.112129,-0.469016,0.0855018,0.0585585,-0.299304,-0.412578,0.43323,
			-0.24062,-0.295829,-0.450792,0.106161,0.0463487,-0.404163,0.136996,-0.0570517,
			-0.433618,-0.125707,-0.250897,0.424875,0.129499,0.378309,0.141674,0.298391,
			-0.064974,0.48114,-0.404042,0.0274824,0.0456458,-0.215657,-0.129197,-0.435307,
			0.0448091,0.336376,-0.354678,-0.32848,-0.431953,0.324012,-0.366029,0.384786,
			0.0147374,0.463636,-0.379505,-0.45171,-0.119848,-0.0872088,-0.0986087,-0.0790029,
			-0.123046,0.407337,0.170162,0.461839,-0.337021,0.248649,-0.125934,-0.0457635,
			-0.461439,0.062432,-0.127688,0.292784,0.295231,-0.117086,-0.24721,-0.157072,
			0.467804,-0.0201916,-0.131672,0.264567,-0.122851,0.400306,-0.316568,-0.131683,
			0.417457,0.0159161,-0.409693,0.235311,-0.495288,0.103123,0.456867,-0.102568,
			0.231551,0.184639,0.478503,-0.296215,0.0933028,0.451563,-0.239728,0.0146776,
			0.136332,-0.0990409,-0.0133887,0.250458,-0.373803,-0.456938,-0.129061,0.193302,
			0.435825,-0.0224244,-0.3709,-0.0161688,0.445597,-0.132256,-0.171514,0.272875,
			-0.202726,-0.322146,0.190799,-0.236055,-0.0422609,0.343692,0.381499,0.20002,
			0.255691,0.474515,-0.0978473,-0.368721,0.224734,0.399518,-0.329299,-0.456971,
			-0.0208408,-0.406064,0.15005,0.452278,-0.0422873,0.0368806,-0.433513,-0.00612992,
			-0.082459,-0.207743,-0.210336,0.253846,-0.403204,-0.423083,0.220915,0.264912,
			0.157945,0.310409,-0.125757,-0.193769,-0.1293,0.206747,-0.331634,0.313721,
			-0.0337717,0.222286,0.494869,-0.137497,0.230828,0.149667,0.181339,-0.492388,
			0.154149,0.445235,0.113271,0.282928,-0.496847,0.296958,0.141817,-0.321521,
			0.0294007,-0.281256,0.0480524,-0.441763,0.0875871,-0.0838964,-0.313555,-0.436089,
			-0.425217,-0.189963,0.444085,0.480727,0.0551155,0.488526,0.19156,-0.258339,
			0.309814,0.434512,-0.371197,0.186827,-0.202755,0.147221,-0.0362412,0.42279,
			-0.258278,0.160157,-0.167713,-0.023247,-0.0312267,0.205929,-0.260051,0.217203,
			0.365197,-0.0890103,-0.0752128,0.454286,0.381419,0.198588,-0.194599,0.329326,
			0.470605,-0.199897,0.498062,-0.0592516,-0.493773,-0.208099,0.191878,-0.00715073,
			-0.416522,-0.304213,0.477595,-0.133801,-0.360568,-0.485244,0.140641,0.237722,
			-0.473243,-0.397847,0.39333,0.263876,0.0826645,0.18543,0.467134,0.0503444,
			0.277197,0.11035,0.488593,-0.451738,0.485389,-0.295264,0.412477,0.165546,
			-0.0376786,-0.451667,-0.0395608,0.299951,-0.210567,0.195102,-0.240721,0.21323,
			0.22037,0.233277,0.122331,0.489779,-0.347622,-0.296677,0.319317,-0.441565,
			0.0385068,-0.309847,0.0994808,-0.207742,-0.408712,0.00676728,0.384132,0.115568,
			-0.453608,0.451924,-0.33097,0.326671,0.111377,0.347299,-0.385877,0.14922,
			-0.385185,-0.0265878,0.183181,-0.366701,-0.0359218,-0.42875,0.0812184,0.0659789,
			-0.244738,-0.261506,-0.484046,-0.115255,0.2573,0.07516,-0.0919145,-0.304292,
			0.0121743,0.213346,0.367386,0.241805,0.494833,0.366748,-0.0141872,0.403305,
			-0.479686,0.136787,0.444056,-0.45922,-0.214328,0.177254,0.0861783,-0.43503};
		int j, iran = i*13;
		double sf = vol->scale[i];
		iran = i*10;
		for(j=0;j<vol->dim[0]*vol->dim[1];j++)
		{
			iran = (iran+1)%397;
#ifdef IGNORE_ZEROS
			if (dat[j])
#endif
				dat[j] += ran[iran]*sf;
		}
	}
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
	*psh = RINT(sh);
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
			mexErrMsgTxt("Arguments must be numeric, real, full and double.");
		}
	}

	vol = get_maps(prhs[0], &n0);
	if (n0!=1)
	{
		free_maps(vol, n0);
		mexErrMsgTxt("Incorrect usage.");
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
	int m0, m1, m2,  n0, n1, n2, nh, sh;
	MAPTYPE *vol, *get_maps();
	double *B0, *B1, *B2, mx, *f, *Alpha, *Beta, *h;

	if (nrhs == 4 && nlhs <=2)
	{
		mexFunction1(nlhs, plhs, nrhs, prhs);
		return;
	}
	if (nrhs != 6 || nlhs > 5) mexErrMsgTxt("Incorrect usage.");

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
