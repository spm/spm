#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

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
void kronutil1(int n1, int n2, int m1, int m2,
	double img[], double b1[], double b2[], double beta[])
{
	int j1,j2, i1,i2;
	double beta1[32];

	/* Zero beta */
	for(j1=0; j1<n1*n2; j1++) beta[j1]=0.0;

	for(i2=0; i2<m2; i2++)
	{
		/* generate small beta */
		for(j1=0; j1<n1; j1++) beta1[j1]=0.0;
		for(i1=0; i1<m1; i1++)
		{
			double wt = img[i1+m1*i2];
			for(j1=0; j1<n1; j1++)
				beta1[j1] += b1[i1+j1*m1]*wt;
		}

		/* kronecker tensor product to increment beta */
		for(j2=0; j2<n2; j2++)
		{
			double wt = b2[i2+j2*m2];
			double *ptrb = beta+(n1*j2);
			for(j1=0; j1<n1; j1++)
				ptrb[j1] += wt*beta1[j1];
		}
	}
}

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
void kronutil2(int n1, int n2, int m1, int m2,
	double img[], double b1[], double b2[], double alpha[])
{
	int j11,j12, j21,j22, i1, i2;
	double alpha1[1024];

	/* Zero alpha */
	for(j21=0; j21<n1*n2; j21++)
		for(j11=0; j11<=j21; j11++)
			alpha[j11+j21*n1*n2]=0.0;

	for(i2=0; i2<m2; i2++)
	{
		/* zero small alpha */
		for(j21=0; j21<n1; j21++)
			for(j11=0; j11<=j21; j11++)
				alpha1[j11+j21*n1]=0.0;

		/* generate upper half of small alpha */
		for(i1=0; i1<m1; i1++)
		{
			double wt = img[i1+m1*i2];
			for(j21=0; j21<n1; j21++)
			{
				double wt2 = wt*b1[i1+j21*m1];
				for(j11=0; j11<=j21; j11++)
					alpha1[j11+j21*n1] += wt2*b1[i1+j11*m1];
			}
		}

		/* kronecker tensor product to increment upper half of large alpha */
		for(j22=0; j22<n2; j22++)
		{
			double wt = b2[i2+j22*m2];
			for(j12=0; j12<=j22; j12++)
			{
				double *ptra = alpha+(n1*j12+n1*n1*n2*j22);
				double wt2 = wt*b2[i2+j12*m2];

				for(j21=0; j21<n1; j21++)
					for(j11=0; j11<=j21; j11++)
						ptra[j11+j21*n1*n2] += wt2*alpha1[j11+j21*n1];
			}
		}
	}
	/* fill in lower triangle of symmetric submatrices of alpha */
	for(j22=0; j22<n2; j22++)
		for(j12=0; j12<j22; j12++)
		{
			double *ptra = alpha+(n1*j12+n1*n1*n2*j22);
			for(j21=0; j21<n1; j21++)
				for(j11=0; j11<j21; j11++)
					ptra[j21+j11*n1*n2] = ptra[j11+j21*n1*n2];
		}

	/* fill in lower triangle of alpha */
	for(j22=0; j22<n2*n1; j22++)
		for(j12=0; j12<j22; j12++)
			alpha[j22+j12*n1*n2]=alpha[j12+j22*n1*n2];
}


void spm_kronutil(int n1, int n2, int m1, int m2,
	double img1[], double img2[], double b1[], double b2[],
	double alpha[], double beta[])
{
	int j11,j12, j21,j22, i1, i2;
	double alpha1[1024], beta1[32];

	/* Zero alpha and beta */
	for(j21=0; j21<n1*n2; j21++)
	{
		beta[j21]=0.0;
		for(j11=0; j11<=j21; j11++)
			alpha[j11+j21*n1*n2]=0.0;
	}

	for(i2=0; i2<m2; i2++)
	{
		/* zero small alpha and beta */
		for(j21=0; j21<n1; j21++)
		{
			beta1[j21]=0.0;
			for(j11=0; j11<=j21; j11++)
				alpha1[j11+j21*n1]=0.0;
		}

		/* generate beta and upper half of small alpha */
		for(i1=0; i1<m1; i1++)
		{
			double pix1 = img1[i1+m1*i2];
			double pix2 = img2[i1+m1*i2];
			for(j21=0; j21<n1; j21++)
			{
				double tmp = pix1*b1[i1+j21*m1];
				beta1[j21] += tmp*pix2;
				for(j11=0; j11<=j21; j11++)
					alpha1[j11+j21*n1] += tmp*pix1*b1[i1+j11*m1];
			}
		}

		/* kronecker tensor product to increment beta and upper half of large alpha */
		for(j22=0; j22<n2; j22++)
		{
			double wt = b2[i2+j22*m2];
			double *ptrb = beta+(n1*j22);

			for(j21=0; j21<n1; j21++)
				ptrb[j21] += wt*beta1[j21];

			for(j12=0; j12<=j22; j12++)
			{
				double *ptra = alpha+(n1*j12+n1*n1*n2*j22);
				double wt2 = wt*b2[i2+j12*m2];

				for(j21=0; j21<n1; j21++)
					for(j11=0; j11<=j21; j11++)
						ptra[j11+j21*n1*n2] += wt2*alpha1[j11+j21*n1];
			}
		}
	}
	/* fill in lower triangle of symmetric submatrices of alpha */
	for(j22=0; j22<n2; j22++)
		for(j12=0; j12<j22; j12++)
		{
			double *ptra = alpha+(n1*j12+n1*n1*n2*j22);
			for(j21=0; j21<n1; j21++)
				for(j11=0; j11<j21; j11++)
					ptra[j21+j11*n1*n2] = ptra[j11+j21*n1*n2];
		}

	/* fill in lower triangle of alpha */
	for(j22=0; j22<n2*n1; j22++)
		for(j12=0; j12<j22; j12++)
			alpha[j22+j12*n1*n2]=alpha[j12+j22*n1*n2];
}


#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int n1, n2, m1, m2;

	double *alpha, *beta, *img, *b1, *b2;

	if (nrhs == 0) mexErrMsgTxt("Incorrect usage");

	if (nrhs != 4) mexErrMsgTxt("4 input arguments required");
	if (nlhs > 1) mexErrMsgTxt("only 1 output argument required");

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("img must be numeric, real, full and double");
	img = mxGetPr(prhs[0]);
	m1 = mxGetM(prhs[0]);
	m2 = mxGetN(prhs[0]);


	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
		mexErrMsgTxt("b1 must be numeric, real, full and double");
	b1 = mxGetPr(prhs[1]);
	if (mxGetM(prhs[1]) != m1)
		mexErrMsgTxt("b1 has incompatible m dimension");
	n1 = mxGetN(prhs[1]);

	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("b2 must be numeric, real, full and double");
	b2 = mxGetPr(prhs[2]);
	if (mxGetM(prhs[2]) != m2)
		mexErrMsgTxt("b2 has incompatible m dimension");
	n2 = mxGetN(prhs[2]);

	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
		mexErrMsgTxt("flag must be numeric, real, full and double");
	if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 1)
		mexErrMsgTxt("flag must contain only one element");

	if (!*mxGetPr(prhs[3]))
	{
		double *beta;
		plhs[0] = mxCreateDoubleMatrix(n1*n2,1, mxREAL);
		beta = mxGetPr(plhs[0]);
		kronutil1(n1,n2,m1,m2,img,b1,b2, beta);
	}
	else
	{
		double *alpha;
		plhs[0] = mxCreateDoubleMatrix(n1*n2,n1*n2, mxREAL);
		alpha = mxGetPr(plhs[0]);
		kronutil2(n1,n2,m1,m2,img,b1,b2, alpha);
	}
}
