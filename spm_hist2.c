#ifndef lint
static char sccsid[] = "%W% John Ashburner %E%";
#endif
#include "math.h"
#define local_rint(a) floor((a)+0.5)

void hist2(double M[16], unsigned char g[], unsigned char f[], const int dg[3], const int df[3], double H[65536], int s[3])
{
	int z;
	for(z=1; z<=dg[2]; z+=s[2])
	{
		int y;
		for(y=1; y<=dg[1]; y+=s[1])
		{
			int x;
			float x3 = M[12] + z*M[ 8] + y*M[4] - 1.0;
			float y3 = M[13] + z*M[ 9] + y*M[5] - 1.0;
			float z3 = M[14] + z*M[10] + y*M[6] - 1.0;
			for(x=1; x<=dg[0]; x+=s[0])
			{
				float x4,y4,z4;
				x4 = x3 + M[0]*x;
				y4 = y3 + M[1]*x;
				z4 = z3 + M[2]*x;
				if (z4>=0 && z4<df[2]-1 && y4>=0 && y4<df[1]-1 && x4>=0 && x4<df[0]-1)
				{
					int k111,k112,k121,k122,k211,k212,k221,k222;
					float dx1, dx2, dy1, dy2, dz1, dz2;
					int off0, ix4, iy4, iz4;
					int vf, vg = g[x-1+dg[0]*(y-1+dg[1]*(z-1))];

					ix4 = floor(x4); dx1=x4-ix4; dx2=1.0-dx1;
					iy4 = floor(y4); dy1=y4-iy4; dy2=1.0-dy1;
					iz4 = floor(z4); dz1=z4-iz4; dz2=1.0-dz1;

					off0 = ix4+df[0]*(iy4+df[1]*iz4);
					k222 = f[off0      ]; k122 = f[off0      +1];
					k212 = f[off0+df[0]]; k112 = f[off0+df[0]+1];
					off0 += df[0]*df[1];
					k221 = f[off0      ]; k121 = f[off0      +1];
					k211 = f[off0+df[0]]; k111 = f[off0+df[0]+1];

					vf = (int)local_rint((((k222*dx2+k122*dx1)*dy2       +
					                 (k212*dx2+k112*dx1)*dy1))*dz2 +
					               (((k221*dx2+k121*dx1)*dy2       +
					                 (k211*dx2+k111*dx1)*dy1))*dz1);

					H[vf+vg*256] ++;

					/*
					H[k222+vg*256] += dx2*dy2*dz2;
					H[k122+vg*256] += dx1*dy2*dz2;
					H[k212+vg*256] += dx2*dy1*dz2;
					H[k112+vg*256] += dx1*dy1*dz2;

					H[k221+vg*256] += dx2*dy2*dz1;
					H[k121+vg*256] += dx1*dy2*dz1;
					H[k211+vg*256] += dx2*dy1*dz1;
					H[k111+vg*256] += dx1*dy1*dz1;
					*/
				}
			}
		}
	}
}

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int *dimsf, *dimsg;
	int s[3];

	if (nrhs>4 || nrhs<3 || nlhs>1) mexErrMsgTxt("Incorrect usage.");

	if (!mxIsNumeric(prhs[0]) || !mxIsUint8(prhs[0]) || mxIsComplex(prhs[0]))
		mexErrMsgTxt("Wrong sort of data (1).");
	if (mxGetNumberOfDimensions(prhs[0]) != 3) mexErrMsgTxt("Wrong number of dims (1).");;
	dimsg  = mxGetDimensions(prhs[0]);

	if (!mxIsNumeric(prhs[1]) || !mxIsUint8(prhs[1]) || mxIsComplex(prhs[1]))
		mexErrMsgTxt("Wrong sort of data (2).");
	if (mxGetNumberOfDimensions(prhs[0]) != 3) mexErrMsgTxt("Wrong number of dims (1).");;
	dimsf  = mxGetDimensions(prhs[1]);

	if (!mxIsNumeric(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
		mexErrMsgTxt("Wrong sort of matrix.");
	if (mxGetM(prhs[2]) != 4 || mxGetN(prhs[2]) != 4)
		mexErrMsgTxt("Matrix must be 4x4.");

	if (nrhs == 4)
	{
		if (!mxIsNumeric(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
		     mxGetM(prhs[3])*mxGetN(prhs[3]) != 3)
			mexErrMsgTxt("Invalid skips.");
		s[0] = mxGetPr(prhs[3])[0];
		s[1] = mxGetPr(prhs[3])[1];
		s[2] = mxGetPr(prhs[3])[2];
	}
	else
	{
		s[0] = s[1] = s[2] = 1;
	}

	plhs[0] = mxCreateDoubleMatrix(256,256,mxREAL);

	hist2(mxGetPr(prhs[2]), (unsigned char *)mxGetPr(prhs[0]), (unsigned char *)mxGetPr(prhs[1]),
		dimsg, dimsf, mxGetPr(plhs[0]), s);

}
