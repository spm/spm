#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif
/*
Generate a residual image
FORMAT VO = spm_resid2(VI,U)
VO     - Output image description
VI     - Vector of input image descriptors
U      - modified design matrix
*/

#include <math.h>
#include "mex.h"
#include "spm_vol_utils.h"

/* A(i,j) = sum_k B(i,k)*C(k,j) */
static void mmult(int ni, int nj, int nk,  double A[], double B[], double C[])
{
	int i, j, k;
	for(i=0; i<ni; i++)
		for(j=0; j<nj; j++)
		{
			A[i+ni*j] = 0.0;
			for(k=0; k<nk; k++)
				A[i+ni*j] += B[i+ni*k]*C[k+nk*j];
		}
}

/* r(k) = sum_i (B(k,i) - sum_j A(i,j)*X(k,j))^2 */
static void resid(int ni, int nj, int nk,  double r[], double A[], double X[], double B[])
{
	int k;
	for(k=0; k<nk; k++)
	{
		double res = 0.0;
		int i, j;
		for(i=0; i<ni; i++)
		{
			double tmp = B[k+i*nk];
			for(j=0; j<nj; j++)
				tmp -= A[i+j*ni]*X[k+j*nk];
			res += tmp*tmp;
		}
		r[k] = res;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *maps, *get_maps();
	double *sptr, *scales, *image, scale, *data1, *data2, *coords;
	short **dptr;
	int j0,j1,j2, nj0,nj1,nj2,nj01,ni,i, maxval, dtype, floatflag;
	double NaN = 0.0/0.0; /* the only way to get a NaN that I know */
	mxArray *wplane_args[3];

	if ((nrhs != 3) || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");

	maps = get_maps(prhs[0], &ni);

	nj0 = maps[0].dim[0];
	nj1 = maps[0].dim[1];
	nj2 = maps[0].dim[2];

	for(i=1; i<ni; i++)
	{
		if (	maps[i].dim[0] != nj0 ||
			maps[i].dim[1] != nj1 ||
			maps[i].dim[2] != nj2)
			{
				free_maps(maps, ni);
				mexErrMsgTxt("Incompatible image dimensions.");
			}
	}

	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) ||
		mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]) ||
		mxGetM(prhs[2]) != ni || mxGetN(prhs[2]) >= ni)
	{
		free_maps(maps, ni);
		mexErrMsgTxt("Weights are badly specified.");
	}


	dtype = get_dtype(prhs[1]);
	if (dtype > 256)
		dtype>>=8;

	if (dtype == 2)
	{
		maxval = 255;
		floatflag = 0;
	}
	else if (dtype == 4)
	{
		maxval = 32767;
		floatflag = 0;
	}
	else if (dtype == 8)
	{
		maxval = 2147483647;
		floatflag = 0;
	}
	else
	{
		floatflag = 1;
	}

	nj2 = nj2;
	nj01 = nj0*nj1;

	/* The compiler doesn't like this line - but I think it's OK */
	wplane_args[0] = prhs[1];
	wplane_args[1] = mxCreateDoubleMatrix(nj0,nj1,mxREAL);
	wplane_args[2] = mxCreateDoubleMatrix(1,1,mxREAL);

	sptr   = mxGetPr(wplane_args[1]);
	image  = (double *)mxCalloc(nj01, sizeof(double));

	data1  = (double *)mxCalloc(nj0*ni, sizeof(double));
	data2  = (double *)mxCalloc(nj0*mxGetN(prhs[2]), sizeof(double));
	coords = (double *)mxCalloc(nj0*3, sizeof(double));


	if (!floatflag)
	{
		scales = (double *)mxCalloc(nj2, sizeof(double));
		dptr   = (short **)mxCalloc(nj2, sizeof(double));
	}

	for(j0=0; j0<nj0; j0++)
		coords[j0+0*nj0] = j0+1.0;

	for(j2=0; j2<nj2; j2++)
	{

		for(j0=0; j0<nj0; j0++)
			coords[j0+2*nj0] = j2+1.0;

		for(j1=0; j1<nj1; j1++)
		{
			for(j0=0; j0<nj0; j0++)
				coords[j0+1*nj0] = j1+1.0;

			for(i=0; i<ni; i++)
				resample(nj0,maps[i],data1+i*nj0,
					coords+0*nj0,coords+1*nj0,coords+2*nj0, 0,NaN);

			mmult(nj0, mxGetN(prhs[2]), ni,  data2, data1, mxGetPr(prhs[2]));

/*
{
	double mx = 0.0;
for(i=0; i<mxGetN(prhs[2])*ni; i++)
	if (mx<data2[i]) mx = data2[i];
printf("data2 - %g\n", mx);
}
*/

			resid(ni, mxGetN(prhs[2]), nj0, 
				sptr+j1*nj0, mxGetPr(prhs[2]), data2, data1);

		}

		if (floatflag)
		{
			mxGetPr(wplane_args[2])[0] = j2+1.0;
			mexCallMATLAB(0, NULL, 3, wplane_args, "spm_write_plane");
		}
		else
		{
			/* Determine maximum and minimum */
			double mx = -9e99;
			for(j1=0; j1<nj01; j1++)
			{
				if (!floatflag && !finite(sptr[j1])) sptr[j1] = 0.0;
				if (sptr[j1]>mx) mx=sptr[j1];
			}

			scales[j2] = mx/32767.0;

			dptr[j2] = (short *)mxCalloc(nj01, sizeof(short));
			for(j1=0; j1<nj01; j1++)
			{
				dptr[j2][j1] = (short)rint(sptr[j1]/scales[j2]);
			}
		}
	}

	if (!floatflag)
	{

		scale = 0.0;
		for(j2=0; j2<nj2; j2++)
		{
			if (scales[j2] > scale)
				scale = scales[j2];
		}
		scale = scale*32767/maxval;

		for(j2=0; j2<nj2; j2++)
		{
			for(j1=0; j1<nj01; j1++)
				sptr[j1] = dptr[j2][j1]*(scales[j2]/scale);

			mxGetPr(wplane_args[2])[0] = j2+1.0;
			mexCallMATLAB(0, NULL, 3, wplane_args, "spm_write_plane");

			mxFree((char *)(dptr[j2]));
		}

		mxFree((char *)scales);
		mxFree((char *)dptr);
	}
	else
	{
		scale = 1.0;
	}

	mxFree((char *)image);

	free_maps(maps, ni);

	if (nlhs == 1)
	{
		plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
		mxGetPr(plhs[0])[0] = scale;
	}
}


static int get_dtype(const mxArray *ptr)
{
	mxArray *tmp;
	double *pr;

	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
		return(NULL);
	}
	tmp=mxGetField(ptr,0,"dim");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find dim.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 4)
	{
		mexErrMsgTxt("Wrong sized dim.");
	}
	pr = mxGetPr(tmp);

	return((int)fabs(pr[3]));
}
