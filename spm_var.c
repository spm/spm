#ifndef lint
static char sccsid[]="%W% %E%";
#endif
 
/*

spm_var.c
% computes to variance from a series of images - a compiled routine
% FORMAT spm_var(N,TYPE,Q,P,[TYPE_O])
% n	-	number of voxels per image
% TYPE	- 	data type (see spm_type.m & volume.h)
% Q	-	filename for variance image
% P	-	matrix of filenames to be averaged (rowwise strings)
% [TYPE_O]	optional type out
%____________________________________________________________________________
%
% spm_var simply computes the variance over a set of images to produce
% an variance image that is written to a named file

*/

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include "mex.h"


#define	N	prhs[0]
#define	TYPE	prhs[1]
#define	Q	prhs[2]
#define	P	prhs[3]
#define	TYPE_O	prhs[4]



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
     double		*n,*b, *bo,*p, k;
     char 		label[512];
     char 		error_msg[512];
     unsigned char      *ku, *ou;
     short 	        *ks, *os;
     int	        *ki, *oi;
     float	        *kf, *of;
     double	        *kd, *od;
     double		*pixel, *pixsq;
     int 		i,j,g,h,check, type, type_out;
     FILE		*fo,*fp[128];

     if (nrhs < 4 || nlhs > 0 || nrhs > 5)
         mexErrMsgTxt("Inappropriate usage.");

     /* open destination file */
     n   = mxGetPr(N); 
     b   = mxGetPr(TYPE); 
     if (nrhs==5) bo  = mxGetPr(TYPE_O); 
     p   = mxGetPr(Q); 
     h   = (int) mxGetN(P);
     g   = (int) mxGetM(P);
     k   = (double) g;

     for (i = 0; i < mxGetN(Q); i++) label[i] = (char) p[i];
     label[i] = (char)0;
     fo  = fopen(label,"w");
     if (fo == (FILE *)0)
         mexErrMsgTxt("Cant open output file.");

     /* open input files */
     p   = mxGetPr(P);

     for (i = 0; i < g; i++) {
	for(j=0;j<512;j++) label[j]= (char)0;
        for (j = 0; j < h; j++) 
	{
          label[j] = (char) p[i + j*g];
	  if (label[j] == (char) 32) label[j] = NULL;
        }


        fp[i]    = fopen(label,"r");
        if ( fp[i] == (FILE *)0) {
	    sprintf(error_msg," \n Cant open %s fd[i] %d i %d ",label,fp[i],i);
	    mexErrMsgTxt(error_msg);
	}

     }

     /* read and accumulate */
     pixel  = (double *) mxCalloc((int) n[0],sizeof(double));
     pixsq =  (double *) mxCalloc((int) n[0],sizeof(double));

     type = (int) b[0];
     if( nrhs == 5)
	     type_out = (int) bo[0];
     else type_out = 64;

     switch(type) 
     {

	case(2) : 

         ku     = (unsigned char *) mxCalloc((int) n[0], sizeof(*ku));
	 if(!ku) mexErrMsgTxt("mxCalloc ku failed ");
         for (j = 0; j < g; j++) {
            check = fread(ku, sizeof(*ku), (int) n[0], fp[j]);
	    if(check != (int) n[0]) mexErrMsgTxt(" check != (int) n[0] ");
            for (i = 0; i < (int) n[0]; i++) {
               pixel[i] += (double) ku[i];
               pixsq[i] += (double) ku[i] * (double) ku[i];
	    }
	    fclose(fp[j]);
         }
	 mxFree(ku);
	 break;


	case(4) : 

         ks     = (short *) mxCalloc((int) n[0], sizeof(*ks));
	 if(!ks) mexErrMsgTxt("mxCalloc ks failed ");
	 
         for (j = 0; j < g; j++) {
            check = fread(ks, sizeof(*ks), (int) n[0], fp[j]);
	    if(check != (int) n[0]) mexErrMsgTxt(" check != (int) n[0] ");
           for (i = 0; i < (int) n[0]; i++) {
               pixel[i] += (double) ks[i];
               pixsq[i] += (double) ks[i] * (double) ks[i];
	    }
	    fclose(fp[j]);
         }
	 mxFree(ks);
	 break;


   	case(8) : 
     
         ki     = (int *) mxCalloc((int) n[0], sizeof(*ki));
	 if(!ki) mexErrMsgTxt("mxCalloc ki failed ");
         for (j = 0; j < g; j++) {
            check = fread(ki, sizeof(*ki), (int) n[0], fp[j]);
	    if(check != (int) n[0]) mexErrMsgTxt(" check != (int) n[0] ");
           for (i = 0; i < (int) n[0]; i++) {
               pixel[i] += (double) ki[i];
               pixsq[i] += (double) ki[i] * (double) ki[i];
	    }
	    fclose(fp[j]);
         }
	 mxFree(ki);
	 break;

   	case(16) : 

         kf     = (float *) mxCalloc((int) n[0], sizeof(*kf));
	 if(!kf) mexErrMsgTxt("mxCalloc kf failed ");
         for (j = 0; j < g; j++) {
            check = fread(kf, sizeof(*kf), (int) n[0], fp[j]);
	    if(check != (int) n[0]) mexErrMsgTxt(" check != (int) n[0] ");
            for (i = 0; i < (int) n[0]; i++) {
               pixel[i] += (double) kf[i];
               pixsq[i] += (double) kf[i] * (double) kf[i];
	    }
	    fclose(fp[j]);
         }
	 mxFree(kf);
	 break;


   	case(64) : 

         kd     = (double *) mxCalloc((int) n[0], sizeof(*kd));
	 if(!kd) mexErrMsgTxt("mxCalloc kd failed ");
         for (j = 0; j < g; j++) {
            check = fread(kd, sizeof(*kd), (int) n[0], fp[j]);
	    if(check != (int) n[0]) mexErrMsgTxt(" check != (int) n[0] ");
            for (i = 0; i < (int) n[0]; i++) {
               pixel[i] += (double) kd[i];
               pixsq[i] += (double) kd[i] * (double) kd[i];
	    }
	    fclose(fp[j]);
         }
	 mxFree(kd);
	 break;


	default : 
             mexErrMsgTxt("Data type not supported");

      } /* switch */
      


     switch(type_out) 
     {

	case(2) : 
      		ou     = (unsigned char *) mxCalloc((int) n[0], sizeof(*ou));
      		if(!ou) mexErrMsgTxt("mxCalloc output failed ");
      		for (i = 0; i < (int) n[0]; i++) 
		ou[i] =  (unsigned char) (pixsq[i]/k - (pixel[i]/k)*(pixel[i]/k));
      		fwrite(ou, sizeof(*ou), (int) n[0], fo);
      		mxFree(ou);

	 break;

	case(4) : 
      		os     = (short *) mxCalloc((int) n[0], sizeof(*os));
      		if(!os) mexErrMsgTxt("mxCalloc output failed ");
      		for (i = 0; i < (int) n[0]; i++) 
           	    os[i] =  (short) (pixsq[i]/k - (pixel[i]/k)*(pixel[i]/k));
      		fwrite(os, sizeof(*os), (int) n[0], fo);
      		mxFree(os);
	 break;

   	case(8) : 
      		oi     = (int *) mxCalloc((int) n[0], sizeof(*oi));
      		if(!oi) mexErrMsgTxt("mxCalloc output failed ");
      		for (i = 0; i < (int) n[0]; i++)
           	    oi[i] =  (int) (pixsq[i]/k - (pixel[i]/k)*(pixel[i]/k));
      		fwrite(oi, sizeof(*oi), (int) n[0], fo);
      		mxFree(oi);
	 break;

   	case(16) : 
      		of     = (float *) mxCalloc((int) n[0], sizeof(*of));
      		if(!of) mexErrMsgTxt("mxCalloc output failed ");
      		for (i = 0; i < (int) n[0]; i++)
           	    of[i] = (float) (pixsq[i]/k - (pixel[i]/k)*(pixel[i]/k));
      		fwrite(of, sizeof(*of), (int) n[0], fo);
      		mxFree(of);
	 break;

   	case(64) : 
      		od     = (double *) mxCalloc((int) n[0], sizeof(*od));
      		if(!od) mexErrMsgTxt("mxCalloc output failed ");
      		for (i = 0; i < (int) n[0]; i++)
           	    od[i] =  (double) (pixsq[i]/k - (pixel[i]/k)*(pixel[i]/k));
      		fwrite(od, sizeof(*od), (int) n[0], fo);
      		mxFree(od);
	 break;


	default : 
             mexErrMsgTxt("Data type not supported");

      } /* switch */
      

      od     = (double *) mxCalloc((int) n[0], sizeof(*od));
      if(!od) mexErrMsgTxt("mxCalloc output failed ");



     /* close files */
     fclose(fo);
}
