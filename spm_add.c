#ifndef lint
static char sccsid[]="%W% Jean-Baptiste Poline, John Ashburner %E%";
#endif
 
/*
spm_mean.c
% averages a series of images - a compiled routine
% FORMAT s = spm_mean(n,TYPE,Q,P,S)
% n	-	number of voxels per image
% TYPE	- 	data type (see spm_type.m & volume.h)
% Q	-	filename for averaged image
% P	-	matrix of filenames to be averaged (rowwise strings)
% S	-	optional vector of scalefactors / weights
% s	-	the scalefactor which should be assigned to Q
%_______________________________________________________________________
%
% spm_mean computes a weighted sum of a set of image files to produce
% an mean image that is written to a named file (Q). No headers are
% read or written, in particular scalefactors are ignored. The image is
% written in the same type as the input images.
%
% The weights (S) default to 1/size(P,1) - resulting an average of the
% image files being written to Q. For a "proper" average including
% scalefactors, the weights s should be specified as sf/length(sf) for
% sf a vector of scalefactors for the image files specified in P.
%
% See also: spm_average - for details on how to usefully use this function
%
%_______________________________________________________________________

*/

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "mex.h"

double rescale(n,d, mx)
int    n;
double d[], mx;
{
	int    i;
	double dmx = 0.0;

	for(i=0; i<n; i++)
		if (d[i] > dmx) dmx = d[i];
	dmx = mx/dmx;
	for(i=0; i<n; i++)
		d[i] *= dmx;
	return(1.0/dmx);
}

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
     double		*p;
     char 		label[1024];
     char 		error_msg[1024];
     unsigned char      *ku;
     short 	        *ks;
     int	        *ki;
     float	        *kf;
     double		*pixel;
     double	        *kd;
     double		scale;
     int 		i,j,g,h, n,b;
     int		fo, fp[1024];

     if ((nrhs != 4 && nrhs != 5) || nlhs > 1)
         mexErrMsgTxt("Inappropriate usage.");

     /* open destination file */
     n   = (int)mxGetPr(prhs[0])[0];
     b   = (int)mxGetPr(prhs[1])[0];
     p   = mxGetPr(prhs[2]);
     h   = (int) mxGetN(prhs[3]);
     g   = (int) mxGetM(prhs[3]);

     if (nrhs == 5)
     {
          if (mxGetN(prhs[4])*mxGetM(prhs[4]) != mxGetM(prhs[3]))
              mexErrMsgTxt("Incorrect number of scalefactors passed.");
     }

     for (i = 0; i < mxGetN(prhs[2]); i++) label[i] = (char) p[i];
     label[i] = 0;
     fo  = open(label, O_RDWR|O_CREAT, 0644);
     if (fo == -1)
         mexErrMsgTxt("Cant open output file.");

     /* open input files */
     p   = mxGetPr(prhs[3]);

     for (i = 0; i < g; i++) {
	for(j=0;j<128;j++) label[j]= 0;
        for (j = 0; j < h; j++){
          label[j] = (char) p[i + j*g];
	/*  mexPrintf("%c%d ",label[j],label[j]); */
	  if (label[j] == (char) 32) label[j] = 0;
        }
	/* mexPrintf(" \n mean :  %s %d ",label,strlen(label)); */

        fp[i]    = open(label,O_RDONLY);
        if (fp[i] == -1) {
	    sprintf(error_msg," \n Cant open %s fd[i] %d i %d ",label,fp[i],i);
	    mexErrMsgTxt(error_msg);
	}

     }

     /* read and accumulate */
     pixel = (double *) mxCalloc((int) (n),sizeof(double));

     if ((int) b == 2) {
         ku     = (unsigned char *) mxCalloc((int) n, sizeof(unsigned char));
         for (j = 0; j < g; j++) {
            if (nrhs == 5) scale = mxGetPr(prhs[4])[j]; else scale = 1.0/((double) g);

            read(fp[j], ku, sizeof(unsigned char)*(int) n);
            for (i = 0; i < (int) n; i++)
               pixel[i] += (double) ku[i]*scale;
	    close(fp[j]);
         }
         scale = rescale(n,pixel, 255.0);
         for (i = 0; i < (int) n; i++)
             ku[i] = (unsigned char) pixel[i];
         write(fo, ku, sizeof(unsigned char)*(int) n);
      }

      else if (b == 4) {
         ks     = (short *) mxCalloc(n, sizeof(short));
         for (j = 0; j < g; j++) {
            if (nrhs == 5) scale = mxGetPr(prhs[4])[j]; else scale = 1.0/((double) g);

            read(fp[j], ks, sizeof(short)*n);
            for (i = 0; i < n; i++)
               pixel[i] += (double) ks[i]*scale;
	    close(fp[j]);
         }
         scale = rescale(n,pixel, 32767.0);
         for (i = 0; i < n; i++)
             ks[i] = (short) pixel[i];
         write(fo, ks, sizeof(short)*n);
      }

      else if (b == 8) {
         ki     = (int *) mxCalloc(n, sizeof(int));
         for (j = 0; j < g; j++) {
            if (nrhs == 5) scale = mxGetPr(prhs[4])[j]; else scale = 1.0/((double) g);

            read(fp[j], ki, sizeof(int)*n);
            for (i = 0; i < n; i++)
               pixel[i] += (double) ki[i]*scale;
	    close(fp[j]);
         }
         scale = rescale(n,pixel, 2147483647.0);
         for (i = 0; i < n; i++)
             ki[i] = (int) pixel[i];
         write(fo, ki, sizeof(int)*n);
      }

      else if (b == 16) {
         kf     = (float *) mxCalloc(n, sizeof(float));
         for (j = 0; j < g; j++) {
            if (nrhs == 5) scale = mxGetPr(prhs[4])[j]; else scale = 1.0/((double) g);

            read(fp[j], kf, sizeof(float)*n);
            for (i = 0; i < n; i++)
               pixel[i] += (double) kf[i]*scale;
	    close(fp[j]);
         }
         scale = 1.0;
         for (i = 0; i < n; i++)
             kf[i] = (float) pixel[i];
         write(fo, kf, sizeof(float)*n);
      }

      else if (b == 32) {
         kd     = (double *) mxCalloc(n, sizeof(double));
         for (j = 0; j < g; j++) {
            if (nrhs == 5) scale = mxGetPr(prhs[4])[j]; else scale = 1.0/((double) g);

            read(fp[j], kd, sizeof(double)*n);
            for (i = 0; i < n; i++)
               pixel[i] += (double) kd[i]*scale;
	    close(fp[j]);
         }
         scale = 1.0;
         for (i = 0; i < n; i++)
             kd[i] = (double) pixel[i];
         write(fo, kd, sizeof(double)*n);
      }

      else  {
             mexErrMsgTxt("Data type not supported");
      }

     close(fo);

     if (nlhs == 1)
     {
          plhs[0] = mxCreateFull(1,1, REAL);
          mxGetPr(plhs[0])[0] = scale;
     }
}
