/*
Matlab MEX file to read an ANT-MPI format continuous file (*.cnt)
this version returns the data in a structure, using fieldnames
similar to the read_ns_hdr and read_asa_msr functions.
*/

/*
Copyright (C) 2003-2004, Robert Oostenveld
Radboud University, Nijmegen, The Netherlands. http://www.kun.nl/fcdonders/
Aalborg University, Denmark. http://www.smi.auc.dk/
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
- Neither the name of the University Nijmegen or the University
  Aalborg, nor the names of its contributors may be used to endorse
  or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include "mex.h"        /* Matlab specific  */
#include "matrix.h"     /* Matlab specific  */
#include "eepmisc.h"    /* MPI-ANT specific */
#include "cnt.h"        /* MPI-ANT specific */

#define NAME "read_eep_cnt"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  /* these variables are MPI-ANT specific for reading the data */
  char filename[256];
  int chan, sample, sb, se, length, chanc, samplec;
  float period;
  int status;
  FILE *fp;
  sraw_t *buf;
  eeg_t *hdr; 

  /* these variables are Matlab specific for interfacing */
  double *ptr;
  const int dims[] = {1, 1};
  mxArray *label;
  mxArray *rate;    
  mxArray *npnt;
  mxArray *nchan;
  mxArray *nsample;
  mxArray *time;
  mxArray *data;

  const int nfields = 7;
  const char *field_names[] = {
    "label",                     /* label           */
    "rate",                      /* 1/(period)          */
    "npnt",                      /* nsample of this segment */
    "nchan",                     /* chanc           */
    "nsample",                   /* nsample of the whole data   */
    "time",                      /*                 */
    "data"};                     /* data            */

  if (nrhs!=3)
    mexErrMsgTxt ("Invalid number of input arguments");

  mxGetString(prhs[0], filename, 256);
  sb = (int)(mxGetScalar(prhs[1]) + 0.5) - 1;
  se = (int)(mxGetScalar(prhs[2]) + 0.5) - 1;
  length = se-sb+1;

  /* open the data file */
 if ((fp = fopen(filename, "rb"))==NULL)
    mexErrMsgTxt ("Could not open file");

  /* read header information */
  hdr = cnt_init_file(filename, fp, &status);
  if (status!=CNTERR_NONE)
    mexErrMsgTxt ("Error reading header from file");

  chanc   = get_cnt_chanc(hdr);
  period  = get_cnt_period(hdr);
  samplec = get_cnt_samplec(hdr);

  if (sb<0)
    mexErrMsgTxt ("Begin sample should be 1 or larger");

  if (se<sb)
    mexErrMsgTxt ("End sample should be similar to, or larger than the begin sample");

  if (chanc<1)
    mexErrMsgTxt ("Invalid number of channels in the data");

  if (se>samplec)
    mexErrMsgTxt ("End sample should be less than the number of samples in the data");

  rate      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(rate     ) = (double)1/period;
  npnt      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(npnt     ) = (double)length;
  nchan     = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(nchan    ) = (double)chanc;
  nsample   = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(nsample  ) = (double)samplec;
  data      = mxCreateDoubleMatrix(chanc,length,mxREAL);
  time      = mxCreateDoubleMatrix(1,length,mxREAL);
  label     = mxCreateCellMatrix(chanc,1);

  ptr = mxGetPr(time);
  for (sample=0; sample<length; sample++)
    ptr[sample] = (double)1000*(sb+sample)*period; 

  for (chan=0; chan<chanc; chan++)
    mxSetCell(label,chan,mxCreateString(get_chan_lab(hdr, chan)));

  /* allocate memory for the data and read it from file */
  buf = (sraw_t *)malloc(CNTBUF_SIZE(hdr, chanc));
  ptr = mxGetPr(data);
  cntseek(hdr, sb); 
  for (sample=0; sample<length; sample++)
    if (cntread(hdr, buf, 1) != CNTERR_NONE)
      mexErrMsgTxt ("Error reading raw data from file");
    else
      for (chan=0; chan<chanc; chan++)
        ptr[sample*chanc+chan] = get_chan_scale(hdr, chan) * buf[chan];

  /* create the struct array with dimensions 1x1 */
  plhs[0] = mxCreateStructArray(2, dims, nfields, field_names);

  /* fill the struct array with the variables */
  mxSetField(plhs[0], 0, "rate",      rate     );
  mxSetField(plhs[0], 0, "npnt",      npnt     );
  mxSetField(plhs[0], 0, "nchan",     nchan    );
  mxSetField(plhs[0], 0, "nsample",   nsample  );
  mxSetField(plhs[0], 0, "data",      data     );
  mxSetField(plhs[0], 0, "time",      time     );
  mxSetField(plhs[0], 0, "label",     label    );

  /* close the file */
  cntclose(hdr);
  fclose(fp);
  free(buf);
 
  return;
}
