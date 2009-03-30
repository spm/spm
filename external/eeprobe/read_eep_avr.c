/*
Matlab MEX file to read an ANT-MPI format average file (*.avr)
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
#include "avr.h"        /* MPI-ANT specific */

#define NAME "read_eep_avr"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  /* these variables are MPI-ANT specific for reading the data */
  char filename[256];
  int chan, sample;
  float *dat;
  float *var;
  double *p1, *p2;
  FILE *fp;
  avr_t hdr;

  /* these variables are Matlab specific for interfacing */
  const int dims[] = {1, 1};
  mxArray *label;
  mxArray *rate;    
  mxArray *npnt;
  mxArray *nchan;
  mxArray *nsweeps;
  mxArray *xmin;
  mxArray *xmax;
  mxArray *time;
  mxArray *data;
  mxArray *variance;
  mxArray *condlab;
  mxArray *condcol;
  mxArray *trialc;
  mxArray *rejtrialc;
  const int nfields = 14;
  const char *field_names[] = {
    "label",                     /* label                   */
    "rate",                      /* 1/(period)              */
    "npnt",                      /* samplec                 */
    "nchan",                     /* chanc                   */
    "nsweeps",                   /* mtrialc                 */
    "xmin",                      /* sample0*period          */
    "xmax",                      /* xmin + samplec*period   */
    "time",                      /*                         */
    "data",                      /* data                    */
    "variance",                  /* variance                */
    "condlab",                   /* condlab                 */
    "condcol",                   /* condcol                 */
    "trialc",                    /* trialc                  */
    "rejtrialc"};                /* rejtrialc               */

  if (nrhs!=1)
    mexErrMsgTxt ("Invalid number of input arguments");
  mxGetString(prhs[0], filename, 256);

  /* open the data file */
  if ((fp = fopen(filename, "rb"))==NULL)
    mexErrMsgTxt ("Could not open file");
 
  /* read header information */
  avropen(&hdr, fp);

  /* allocate memory for the average data and read it from file */
  dat = (float *)malloc(hdr.chanc * hdr.samplec * sizeof(float));
  for (chan=0; chan<hdr.chanc; chan++)
  {
    avrseek(&hdr, fp, chan, AVRBAND_MEAN);
    avrread(fp, dat+chan*hdr.samplec, hdr.samplec);
  }

  /* allocate memory for the variance data and read it from file */
  var = (float *)malloc(hdr.chanc * hdr.samplec * sizeof(float));
  for (chan=0; chan<hdr.chanc; chan++)
  {
    avrseek(&hdr, fp, chan, AVRBAND_VAR);
    avrread(fp, var+chan*hdr.samplec, hdr.samplec);
  }

  /* create matlab variables and assign the numerical values of the header */
  rate      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(rate     ) = (double)1/hdr.period;
  npnt      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(npnt     ) = (double)hdr.samplec;
  nchan     = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(nchan    ) = (double)hdr.chanc;
  nsweeps   = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(nsweeps  ) = (double)hdr.mtrialc;
  xmin      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(xmin     ) = (double)1000*hdr.sample0*hdr.period;
  xmax      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(xmax     ) = (double)1000*(hdr.sample0*hdr.period + hdr.samplec*hdr.period);
  trialc    = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(trialc   ) = (double)hdr.trialc;
  rejtrialc = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(rejtrialc) = (double)hdr.rejtrialc;
  condlab   = mxCreateString(hdr.condlab);
  condcol   = mxCreateString(hdr.condcol);
  data      = mxCreateDoubleMatrix(hdr.chanc,hdr.samplec,mxREAL);
  variance  = mxCreateDoubleMatrix(hdr.chanc,hdr.samplec,mxREAL);
  time      = mxCreateDoubleMatrix(1,hdr.samplec,mxREAL);
  label     = mxCreateCellMatrix(hdr.chanc,1);

  /* create matlab variables and assign the data and variance for each sample */
  p1 = mxGetPr(data);
  p2 = mxGetPr(variance);

  for (chan=0; chan<hdr.chanc; chan++)
  {
    mxSetCell(label,chan,mxCreateString(hdr.chanv[chan].lab));
    for (sample=0; sample<hdr.samplec; sample++)
    {
      *(mxGetPr(time)+sample) = (double)1000*(hdr.sample0+sample)*hdr.period;
      p1[sample*hdr.chanc+chan] = dat[chan*hdr.samplec+sample];
      p2[sample*hdr.chanc+chan] = var[chan*hdr.samplec+sample];
    }
  }

  /* create the struct array with dimensions 1x1 */
  plhs[0] = mxCreateStructArray(2, dims, nfields, field_names);

  /* fill the struct array with the variables */
  mxSetField(plhs[0], 0, "rate",      rate     );
  mxSetField(plhs[0], 0, "npnt",      npnt     );
  mxSetField(plhs[0], 0, "nchan",     nchan    );
  mxSetField(plhs[0], 0, "nsweeps",   nsweeps  );
  mxSetField(plhs[0], 0, "xmin",      xmin     );
  mxSetField(plhs[0], 0, "xmax",      xmax     );
  mxSetField(plhs[0], 0, "trialc",    trialc   );
  mxSetField(plhs[0], 0, "rejtrialc", rejtrialc);
  mxSetField(plhs[0], 0, "condlab",   condlab  );
  mxSetField(plhs[0], 0, "condcol",   condcol  );
  mxSetField(plhs[0], 0, "data",      data     );
  mxSetField(plhs[0], 0, "variance",  variance );
  mxSetField(plhs[0], 0, "time",      time     );
  mxSetField(plhs[0], 0, "label",     label    );

  /* close the data file */
  avrclose(&hdr);
  fclose(fp);

  return;
}
