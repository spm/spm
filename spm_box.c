#ifndef lint
static char sccsid[]="%W% %E%";
#endif
/*

spm_box.c
% integrates a volume image iver x, y and z - a compiled routine
% FORMAT [X Y Z] = spm_box(P,DIM,TYPE)
% P	-	filename (unsigned char)
% DIM	-	[x y z] - image size {voxels}
% TYPE  -       data type (see spm_type.m & volume.h)
% X,Y,Z	-	integrated 1-dimensional images
%____________________________________________________________________________
%
% spm_box simply integrates (sums) all voxel values over x, y and z to
% give three vectors of integrated voxel values.  It is used primarily
% determine the bounding box in which the 'object' lies in image space
%
% see also spm_bb.m

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

#define SIZE    sizeof(char)

#ifdef __STDC__
void mexFunction(
	int		nlhs,
	Matrix	*plhs[],
	int		nrhs,
	Matrix	*prhs[]
	)
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
     double		*n,*p,*d,*X,*Y,*Z,*t;
     char 		label[128];
     unsigned char      c;
     int 		i,j,x,y,z,N;
     FILE		*f;
     unsigned char      *ku;
     short 	        *ks;
     int	        *ki;
     float	        *kf;
     double	        *kd;


     if (nrhs != 3)
         mexErrMsgTxt("Inappropriate usage.");

     /* open destination file */
     p   = mxGetPr(prhs[0]);
     d   = mxGetPr(prhs[1]);
     t   = mxGetPr(prhs[2]);
     N   = (int) d[0]*d[1]*d[2];

     for (i = 0; i < mxGetN(prhs[0]); i++) label[i] = (char) p[i];
     label[i] = NULL;
     f  = fopen(label,"r");
     if (f == (FILE *)0)
            mexErrMsgTxt("Cant open input file.");


     /* create output matrices */
     plhs[0] = mxCreateFull((int) 1,(int) d[0],REAL);
     plhs[1] = mxCreateFull((int) 1,(int) d[1],REAL);
     plhs[2] = mxCreateFull((int) 1,(int) d[2],REAL);
     X       = mxGetPr(plhs[0]);
     Y       = mxGetPr(plhs[1]);
     Z       = mxGetPr(plhs[2]);
   

     /* integrate */
     if ((int) t[0] == 2) {
         ku     = (unsigned char *) mxCalloc(N, sizeof(unsigned char));
         fread(ku, sizeof(unsigned char), N, f);
         for (i = 0; i < N; i++) {
             z     = (int) floor(i / (d[0]*d[1]));
             y     = (int) floor((i - z*(d[0]*d[1]) ) / (d[0]));
             x     = (int) i - z*(d[0]*d[1]) - y*d[0];
             X[x] += (double) ku[i];
             Y[y] += (double) ku[i];
             Z[z] += (double) ku[i];
	 }
     }
     else if ((int) t[0] == 4) {
         ks     = (short *) mxCalloc(N, sizeof(short));
         fread(ks, sizeof(short), N, f);
         for (i = 0; i < N; i++) {
             z     = (int) floor(i / (d[0]*d[1]));
             y     = (int) floor((i - z*(d[0]*d[1]) ) / (d[0]));
             x     = (int) i - z*(d[0]*d[1]) - y*d[0];
             X[x] += (double) ks[i];
             Y[y] += (double) ks[i];
             Z[z] += (double) ks[i];
	 }
     }
     else if ((int) t[0] == 8) {
         ki     = (int *) mxCalloc(N, sizeof(int));
         fread(ki, sizeof(int), N, f);
         for (i = 0; i < N; i++) {
             z     = (int) floor(i / (d[0]*d[1]));
             y     = (int) floor((i - z*(d[0]*d[1]) ) / (d[0]));
             x     = (int) i - z*(d[0]*d[1]) - y*d[0];
             X[x] += (double) ki[i];
             Y[y] += (double) ki[i];
             Z[z] += (double) ki[i];
	 }
     }
     else if ((int) t[0] == 16) {
         kf     = (float *) mxCalloc(N, sizeof(float));
         fread(kf, sizeof(float), N, f);
         for (i = 0; i < N; i++) {
             z     = (int) floor(i / (d[0]*d[1]));
             y     = (int) floor((i - z*(d[0]*d[1]) ) / (d[0]));
             x     = (int) i - z*(d[0]*d[1]) - y*d[0];
             X[x] += (double) kf[i];
             Y[y] += (double) kf[i];
             Z[z] += (double) kf[i];
	 }
     }
     else if ((int) t[0] == 64) {
         kd     = (double *) mxCalloc(N, sizeof(double));
         fread(kd, sizeof(double), N, f);
         for (i = 0; i < N; i++) {
             z     = (int) floor(i / (d[0]*d[1]));
             y     = (int) floor((i - z*(d[0]*d[1]) ) / (d[0]));
             x     = (int) i - z*(d[0]*d[1]) - y*d[0];
             X[x] += (double) kd[i];
             Y[y] += (double) kd[i];
             Z[z] += (double) kd[i];
	 }
     }
     else
	  mexErrMsgTxt("Data type not supported.");
  
     /* close files */
     fclose(f);
}
