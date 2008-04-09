#include <math.h>
#include "mex.h"
#include "matrix.h"

/*
 * Copyright (C) 2006, Robert Oostenveld, F.C. Donders Ccentre for Cognitive Neuroimaging
 * Copyright (C) 2007, Peter Desain, Nijmegen Institute for Cognition and Information
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * $Log: read_24bit.c,v $
 * Revision 1.4  2008/04/09 10:10:20  roboos
 * corrected the text of an error message
 *
 * Revision 1.3  2007/10/02 08:49:37  roboos
 * added brackets around the bitshift operations
 *
 * Revision 1.2  2007/09/12 12:43:39  roboos
 * added peter to copyright, merged suggestions from nici (see email 20070821), changed to using bitshift, changed the handing of the sign bit
 *
 *
 */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  FILE *fp;
  mxArray *dat;
  double *dat_p;
  char *filename, *buf;
  long int offset, numwords, fnlen, i, j, num, b1, b2, b3;
  
  if (nrhs != 3)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  
  /* get the filename */
  fnlen    = mxGetN      (prhs[0]);
  filename = mxCalloc(fnlen+1, sizeof(char));
  mxGetString(prhs[0], filename, fnlen+1);
  
  /* get values and typecast to integer */
  offset   = (long int)(mxGetScalar (prhs[1]));
  numwords = (long int)(mxGetScalar (prhs[2]));

  /*
  printf("filename = %s\n", filename);
  printf("fnlen    = %d\n", fnlen);
  printf("offset   = %d\n", offset);
  printf("numwords = %d\n", numwords);
  */
  
  /* open the file and read the desired bytes */
  fp = fopen(filename, "rb");
  if (fp==NULL)
  {
    printf("error opening file: %s\n", filename);
    return;
  }
  fseek(fp, offset, SEEK_SET);

  buf = mxCalloc(3*numwords, sizeof(char));
  num = fread(buf, sizeof(char), 3*numwords, fp);
  if (num<3*numwords)
  {
    printf("error reading from %s\n", filename);
    fclose(fp);
    return;
  }
  else
  {
    fclose(fp);
  }

  /* convert every thee bytes into one word */
  dat   = mxCreateDoubleMatrix (1, numwords, mxREAL);
  dat_p = mxGetData (dat);

  for (i=0; i<numwords; i++)
  {
    j = i*3;
    b1 = (0x000000FF & ((long int)buf[j]));
    b2 = (0x000000FF & ((long int)buf[j+1]));
    b3 = (0x000000FF & ((long int)buf[j+2]));
    dat_p[i] = ((long int) ((b3 << 24) | (b2 << 16) | (b1 << 8)))/256;
  }

  /* assign the output parameters */
  plhs[0] = dat;
  return;
}

