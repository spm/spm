#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include "cmex.h"

/* C = A'*A */
void atranspa(m,n,A,C)
int m,n;
double A[/* m,n */], C[/* n,n */];
{
	int i, j1,j2;
	double *p1, *p2, c;

	/* Generate half of symmetric matrix C */
	for (j1=0;j1<n;j1++)
	{
		p1 = &(A[j1*m]);
		for (j2=0;j2<=j1;j2++)
		{
			p2 = &(A[j2*m]);
			c = 0.0;

			/* Work down columns in inner loop
			   to reduce paging */
			for(i=0; i<m; i++)
				c += p1[i]*p2[i];

			C[j1*n+j2] = c;
		}
	}

	/* Generate other half */
	for(j1=0; j1<n; j1++)
		for (j2=0;j2<j1;j2++)
			C[j2*n+j1] = C[j1*n+j2];
}


#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	unsigned int n;
	unsigned int m;
	double *C;
	double *A;

	if (nrhs == 0) mexErrMsgTxt("Incorrect usage.");
	if (nrhs != 1) mexErrMsgTxt("Only 1 input argument required.");
	if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required.");

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsFull(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("spm_atranspa: A must be numeric, real, full and double");
	A = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	plhs[0] = mxCreateFull(n,n, REAL);
	C = mxGetPr(plhs[0]);

	atranspa(m,n,A,C);

}

