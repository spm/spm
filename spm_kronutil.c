/*
Performs:

	[m1,m2]=size(img1);
	[m1,m2]=size(img2);
	[m1,n1]=size(b1);
	[m2,n2]=size(b2);

	alpha = zeros(n1*n2,n1*n2);
	beta  = zeros(n1*n2,1);
	for i=1:m2
		tmp   = kron(img1(:,i),ones(1,n1)).*b1;
		alpha = alpha + kron(b2(i,:)'*b1(i,:),  tmp'*tmp);
		beta  = beta  + kron(b2(i,:)', tmp'*img2(:,i));
	end

which is equivalent to, but a lot faster than:

	B     = kron(b2,b1);
	A     = diag(img1(:))*B;
	b     = img2(:);
	alpha = A'*A;
	beta  = A'*b;

*/

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


#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int n1;
	unsigned int m2;
	unsigned int n2;
	unsigned int m1;
	double *alpha;
	double *beta;
	double *img1;
	double *img2;
	double *b1;
	double *b2;

	if (nrhs == 0) mexErrMsgTxt("usage: [alpha,beta]=spm_kronutil(img1,img2,b1,b2)");
	if (nrhs != 4) mexErrMsgTxt("spm_kronutil: 4 input arguments required");
	if (nlhs > 2) mexErrMsgTxt("spm_kronutil: only 2 output arguments required");

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("spm_kronutil: img1 must be numeric, real, full and double");
	img1 = mxGetPr(prhs[0]);
	m1 = mxGetM(prhs[0]);
	m2 = mxGetN(prhs[0]);
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
		mexErrMsgTxt("spm_kronutil: img2 must be numeric, real, full and double");
	img2 = mxGetPr(prhs[1]);
	if (mxGetM(prhs[1]) != m1)
		mexErrMsgTxt("spm_kronutil: img2 has incompatible m dimension");
	if (mxGetN(prhs[1]) != m2)
		mexErrMsgTxt("spm_kronutil: img2 has incompatible no of columns");
	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("spm_kronutil: b1 must be numeric, real, full and double");
	b1 = mxGetPr(prhs[2]);
	if (mxGetM(prhs[2]) != m1)
		mexErrMsgTxt("spm_kronutil: b1 has incompatible m dimension");
	n1 = mxGetN(prhs[2]);
	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
		mexErrMsgTxt("spm_kronutil: b2 must be numeric, real, full and double");
	b2 = mxGetPr(prhs[3]);
	if (mxGetM(prhs[3]) != m2)
		mexErrMsgTxt("spm_kronutil: b2 has incompatible m dimension");
	n2 = mxGetN(prhs[3]);

	plhs[0] = mxCreateDoubleMatrix(n1*n2,n1*n2, mxREAL);
	alpha = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(n1*n2,1, mxREAL);
	beta = mxGetPr(plhs[1]);

	(void)spm_kronutil(n1,n2,m1,m2,img1,img2,b1,b2,alpha,beta);

}
