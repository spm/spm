#ifndef lint
static char sccsid[]="%W% %E%";
#endif
 
/*

spm_hist.c
% returns the histogram of volxe values - a compiled routine
% FORMAT [Y X] = spm_hist(P,DIM,TYPE)
% P	-	filename (unsigned char)
% DIM	-	[x y z] - image size {voxels}
% TYPE  -       data type (see spm_type.m & volume.h)
% X	-	Bins
% Y     -       Frequency
%____________________________________________________________________________
%
% spm_hist simply returns the histogram of volxe values in image (P)
% over 64 bins

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

#define SIZE    sizeof(char)
#define BINS    64

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
     double		*n,*p,*d,*X,*Y,*t;
     double		mx;
     char 		label[128];
     unsigned char      c;
     int 		i,j,N;
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
     plhs[0] = mxCreateFull((int) 1,(int) BINS,REAL);
     plhs[1] = mxCreateFull((int) 1,(int) BINS,REAL);
     Y       = mxGetPr(plhs[0]);
     X       = mxGetPr(plhs[1]);
   

     /* integrate */
     if ((int) t[0] == 2) {
         ku     = (unsigned char *) mxCalloc(N, sizeof(unsigned char));
         fread(ku, sizeof(unsigned char), N, f);
	 mx = (double) ku[0];
         for (i = 0; i < N; i++)
	 if (mx < (double) ku[i]) mx = (double) ku[i];
         for (i = 0; i < N; i++) {
		j  = (int) ((BINS - 1)*((double) ku[i])/mx);
		if (j >= 0) Y[j]++ ; 
 	 }
     }
     else if ((int) t[0] == 4) {
         ks     = (short *) mxCalloc(N, sizeof(short));
         fread(ks, sizeof(short), N, f);
	 mx = (double) ks[0];
         for (i = 0; i < N; i++)
	 if (mx < (double) ks[i]) mx = (double) ks[i];
         for (i = 0; i < N; i++) {
		j  = (int) ((BINS - 1)*((double) ks[i])/mx);
		if (j >= 0) Y[j]++ ; 
	 }
     }
     else if ((int) t[0] == 8) {
         ki     = (int *) mxCalloc(N, sizeof(int));
         fread(ki, sizeof(int), N, f);
	 mx = (double) ki[0];
         for (i = 0; i < N; i++)
	 if (mx < (double) ki[i]) mx = (double) ki[i];
         for (i = 0; i < N; i++) {
		j  = (int) ((BINS - 1)*((double) ki[i])/mx);
		if (j >= 0) Y[j]++ ; 
	 }
     }
     else if ((int) t[0] == 16) {
         kf     = (float *) mxCalloc(N, sizeof(float));
         fread(kf, sizeof(float), N, f);
	 mx = (double) kf[0];
         for (i = 0; i < N; i++)
	 if (mx < (double) kf[i]) mx = (double) kf[i];
         for (i = 0; i < N; i++) {
		j  = (int) ((BINS - 1)*((double) kf[i])/mx);
		if (j >= 0) Y[j]++ ; 
	 }
     }
     else if ((int) t[0] == 64) {
         kd     = (double *) mxCalloc(N, sizeof(double));
         fread(kd, sizeof(double), N, f);
  	 mx = (double) kd[0];
         for (i = 0; i < N; i++)
	 if (mx < (double) kd[i]) mx = (double) kd[i];
         for (i = 0; i < N; i++) {
		j  = (int) ((BINS - 1)*((double) kd[i])/mx);
		if (j >= 0) Y[j]++ ; 
	 }
     }
     else
	  mexErrMsgTxt("Data type not supported.");
  
     /* fill in bins */
     for (i = 0; i < BINS; i++) X[i] = (double) i * mx/(BINS - 1);

     /* close files */
     fclose(f);
}
