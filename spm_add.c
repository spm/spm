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
     char 		label[128];
     char 		error_msg[128];
     unsigned char      *ku;
     short 	        *ks;
     int	        *ki;
     float	        *kf;
     double		*pixel;
     double	        *kd;
     int 		i,j,g,h;
     FILE		*fo,*fp[128];

     if (nrhs != 4 || nlhs > 0)
         mexErrMsgTxt("Inappropriate usage.");

     /* open destination file */
     n   = mxGetPr(prhs[0]);
     b   = mxGetPr(prhs[1]);
     p   = mxGetPr(prhs[2]);
     h   = (int) mxGetN(prhs[3]);
     g   = (int) mxGetM(prhs[3]);

     for (i = 0; i < mxGetN(prhs[2]); i++) label[i] = (char) p[i];
     label[i] = NULL;
     fo  = fopen(label,"w");
     if (fo == (FILE *)0)
         mexErrMsgTxt("Cant open output file.");

     /* open input files */
     p   = mxGetPr(prhs[3]);

     for (i = 0; i < g; i++) {
	for(j=0;j<128;j++) label[j]= NULL;
        for (j = 0; j < h; j++){
          label[j] = (char) p[i + j*g];
	/*  mexPrintf("%c%d ",label[j],label[j]); */
	  if (label[j] == (char) 32) label[j] = NULL;
        }
	/* mexPrintf(" \n mean :  %s %d ",label,strlen(label)); */

        fp[i]    = fopen(label,"r");
        if (fp[i] == (FILE *)0) {
	    sprintf(error_msg," \n Cant open %s fd[i] %d i %d ",label,fp[i],i);
	    mexErrMsgTxt(error_msg);
	}

     }

     /* read and accumulate */
     pixel = (double *) mxCalloc((int) (n[0]),sizeof(double));

     if ((int) b[0] == 2) {
         ku     = (unsigned char *) mxCalloc((int) n[0], sizeof(unsigned char));
         for (j = 0; j < g; j++) {
            fread(ku, sizeof(unsigned char), (int) n[0], fp[j]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) ku[i];
	    fclose(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             ku[i] = (unsigned char) (pixel[i]/((double) g));
         fwrite(ku, sizeof(unsigned char), (int) n[0], fo);
      }

      else if ((int) b[0] == 4) {
         ks     = (short *) mxCalloc((int) n[0], sizeof(short));
         for (j = 0; j < g; j++) {
            fread(ks, sizeof(short), (int) n[0], fp[j]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) ks[i];
	    fclose(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             ks[i] = (short) (pixel[i]/((double) g));
         fwrite(ks, sizeof(short), (int) n[0], fo);
      }

      else if ((int) b[0] == 8) {
         ki     = (int *) mxCalloc((int) n[0], sizeof(int));
         for (j = 0; j < g; j++) {
            fread(ki, sizeof(int), (int) n[0], fp[j]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) ki[i];
	    fclose(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             ki[i] = (int) (pixel[i]/((double) g));
         fwrite(ki, sizeof(int), (int) n[0], fo);
      }

      else if ((int) b[0] == 16) {
         kf     = (float *) mxCalloc((int) n[0], sizeof(float));
         for (j = 0; j < g; j++) {
            fread(kf, sizeof(float), (int) n[0], fp[j]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) kf[i];
	    fclose(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             kf[i] = (float) (pixel[i]/((double) g));
         fwrite(kf, sizeof(float), (int) n[0], fo);
      }

      else if ((int) b[0] == 32) {
         kd     = (double *) mxCalloc((int) n[0], sizeof(double));
         for (j = 0; j < g; j++) {
            fread(kd, sizeof(double), (int) n[0], fp[j]);
            for (i = 0; i < (int) n[0]; i++)
               pixel[i] += (double) kd[i];
	    fclose(fp[j]);
         }
         for (i = 0; i < (int) n[0]; i++)
             kd[i] = (double) (pixel[i]/((double) g));
         fwrite(kd, sizeof(double), (int) n[0], fo);
      }

      else  {
             mexErrMsgTxt("Data type not supported");
      }

     /* close files */
     fclose(fo);
}
