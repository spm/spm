#ifndef lint
static char sccsid[]="%W% JB Poline %E%";
#endif
 
/*

spm_mean.c
% averages a series of images - a compiled routine
% FORMAT spm_mean(n,TYPE,Q,P)
% n	-	number of voxels per image
% TYPE	- 	data type (see spm_type.m & volume.h)
% Q	-	filename for averaged image
% P	-	matrix of filenames to be averaged (rowwise strings)
%____________________________________________________________________________
%
% spm_mean simply averages a set of images to produce an average that
% is written to a named file

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
     double		*n,*b,*p;
     char 		label[1024];
     char 		error_msg[1024];
     unsigned char      *ku;
     short 	        *ks;
     int	        *ki;
     float	        *kf;
     double		*pixel;
     double	        *kd;
     int 		i,j,g,h;
     int		fo, fp[1024];

     if (nrhs != 4 || nlhs > 0)
         mexErrMsgTxt("Inappropriate usage.");

     /* open destination file */
     n   = mxGetPr(prhs[0]);
     b   = mxGetPr(prhs[1]);
     p   = mxGetPr(prhs[2]);
     h   = (int) mxGetN(prhs[3]);
     g   = (int) mxGetM(prhs[3]);

     for (i = 0; i < mxGetN(prhs[2]); i++) label[i] = (char) p[i];
     label[i] = 0;
     fo  = open(label, O_WRONLY|O_CREAT);
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
     pixel = (double *) mxCalloc((int) (n[0]),sizeof(double));

     if ((int) b[0] == 2) {
         ku     = (unsigned char *) mxCalloc((int) n[0], sizeof(unsigned char));
         for (j = 0; j < g; j++) {
            read(fp[j], ku, sizeof(unsigned char)*(int) n[0]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) ku[i];
	    close(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             ku[i] = (unsigned char) (pixel[i]/((double) g));
         write(fo, ku, sizeof(unsigned char)*(int) n[0]);
      }

      else if ((int) b[0] == 4) {
         ks     = (short *) mxCalloc((int) n[0], sizeof(short));
         for (j = 0; j < g; j++) {
            read(fp[j], ks, sizeof(short)*(int) n[0]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) ks[i];
	    close(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             ks[i] = (short) (pixel[i]/((double) g));
         write(fo, ks, sizeof(short)*(int) n[0]);
      }

      else if ((int) b[0] == 8) {
         ki     = (int *) mxCalloc((int) n[0], sizeof(int));
         for (j = 0; j < g; j++) {
            read(fp[j], ki, sizeof(int)*(int) n[0]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) ki[i];
	    close(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             ki[i] = (int) (pixel[i]/((double) g));
         write(fo, ki, sizeof(int)*(int) n[0]);
      }

      else if ((int) b[0] == 16) {
         kf     = (float *) mxCalloc((int) n[0], sizeof(float));
         for (j = 0; j < g; j++) {
            read(fp[j], kf, sizeof(float)*(int) n[0]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) kf[i];
	    close(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             kf[i] = (float) (pixel[i]/((double) g));
         write(fo, kf, sizeof(float)*(int) n[0]);
      }

      else if ((int) b[0] == 32) {
         kd     = (double *) mxCalloc((int) n[0], sizeof(double));
         for (j = 0; j < g; j++) {
            read(fp[j], kd, sizeof(double)*(int) n[0]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) kd[i];
	    close(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             kd[i] = (double) (pixel[i]/((double) g));
         write(fo, kd, sizeof(double)*(int) n[0]);
      }

      else  {
             mexErrMsgTxt("Data type not supported");
      }

     /* close files */
     close(fo);
}
