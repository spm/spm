/********************************************************************************
 *                                                                              *
 * this file is part of:                                                        *
 * libeep, the project for reading and writing avr/cnt eeg and related files    *
 *                                                                              *
 ********************************************************************************
 *                                                                              *
 * LICENSE:Copyright (c) 2003-2009,                                             *
 * Advanced Neuro Technology (ANT) B.V., Enschede, The Netherlands              *
 * Max-Planck Institute for Human Cognitive & Brain Sciences, Leipzig, Germany  *
 *                                                                              *
 ********************************************************************************
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify         *
 * it under the terms of the GNU Lesser General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License, or            *
 * (at your option) any later version.                                          *
 *                                                                              *
 * This library is distributed WITHOUT ANY WARRANTY; even the implied warranty  *
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 * GNU Lesser General Public License for more details.                          *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this program. If not, see <http://www.gnu.org/licenses/>          *
 *                                                                              *
 *******************************************************************************/

#include <math.h>
#include "mex.h"		/* Matlab specific  */
#include "matrix.h"		/* Matlab specific  */
#include "cnt/cnt.h"

#define NAME "read_eep_avr"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int status;
  int has_variance;
  /* these variables are MPI-ANT specific for reading the data */
  char filename[256];
  int chan, sample;
  float *dat;
  float *var;
  double *p1, *p2;
  FILE *fp;
  eeg_t* hdr;

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
  mxArray *psi;
  mxArray *trialc;
  mxArray *rejtrialc;
  const int nfields = 14;
  const char *field_names[] = {
    "label",                     /* label */
    "rate",                      /* 1/(period) */
    "npnt",                      /* samplec */
    "nchan",                     /* chanc */
    "nsweeps",                   /* mtrialc */
    "xmin",                      /* sample0*period */
    "xmax",                      /* xmin + samplec*period */
    "time",                      /*  */
    "data",                      /* data */
    "condlab",                   /* condlab */
    "condcol",                   /* condcol */
    "psi",                       /* pre-stimulus interval */
    "trialc",                    /* trialc */
    "rejtrialc"};                /* rejtrialc */

  if (nrhs!=1)
    mexErrMsgTxt ("Invalid number of input arguments");
  mxGetString(prhs[0], filename, 256);

  /* open the data file */
  if ((fp = eepio_fopen(filename, "rb"))==NULL)
    mexErrMsgTxt ("Could not open file");

  /* read header information */
  /*  avropen(&hdr, fp); */
  hdr = eep_init_from_file(filename, fp, &status);

  has_variance = eep_has_data_of_type(hdr, DATATYPE_STDDEV);

  if( hdr == 0 || status != CNTERR_NONE)
    mexErrMsgTxt ("Could not initialize file");

  /* allocate memory for the average data and read it from file */

  dat = (float*) malloc(FLOAT_CNTBUF_SIZE(hdr, 1));
  if(has_variance)
    var = (float *)malloc(FLOAT_CNTBUF_SIZE(hdr, 1));

  rate      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(rate     ) = (double)eep_get_rate(hdr);
  npnt      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(npnt     ) = (double)eep_get_samplec(hdr);
  nchan     = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(nchan    ) = (double)eep_get_chanc(hdr);
  nsweeps   = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(nsweeps  ) = (double)eep_get_averaged_trials(hdr);
  xmin      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(xmin     ) = (double)1000*eep_get_sample0(hdr)*eep_get_period(hdr);
  xmax      = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(xmax     ) = (double)1000*(eep_get_sample0(hdr)*eep_get_period(hdr) + eep_get_samplec(hdr)*eep_get_period(hdr));
  trialc    = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(trialc   ) = (double)eep_get_total_trials(hdr);
  rejtrialc = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(rejtrialc) = (double)eep_get_total_trials(hdr) - eep_get_averaged_trials(hdr);
  condlab   = mxCreateString(eep_get_conditionlabel(hdr));
  condcol   = mxCreateString(eep_get_conditioncolor(hdr));
  psi       = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(psi      ) = (double)eep_get_pre_stimulus_interval(hdr);
  data      = mxCreateDoubleMatrix(eep_get_chanc(hdr),eep_get_samplec(hdr),mxREAL);
  if(has_variance)
    variance  = mxCreateDoubleMatrix(eep_get_chanc(hdr),eep_get_samplec(hdr),mxREAL);
  time      = mxCreateDoubleMatrix(1,eep_get_samplec(hdr),mxREAL);
  label     = mxCreateCellMatrix(eep_get_chanc(hdr),1);

  p1 = mxGetPr(data);
  if(has_variance)
    p2 = mxGetPr(variance);

  for (chan=0; chan<eep_get_chanc(hdr); chan++)
  {
    mxSetCell(label,chan,mxCreateString(eep_get_chan_label(hdr, chan)));
  }


  for (sample=0; sample<eep_get_samplec(hdr); sample++)
  {
    *(mxGetPr(time)+sample) = (double)1000*(eep_get_sample0(hdr)+sample)*eep_get_period(hdr);
    eep_read_float(hdr, DATATYPE_AVERAGE, dat,1);
    if(has_variance)
      eep_read_float(hdr, DATATYPE_STDDEV, var,1);

    for(chan = 0; chan < eep_get_chanc(hdr); chan++)
    {
      p1[sample*eep_get_chanc(hdr)+chan] = dat[chan];
      if(has_variance)
        p2[sample*eep_get_chanc(hdr)+chan] = var[chan];
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
  mxSetField(plhs[0], 0, "psi",       psi      );
  mxSetField(plhs[0], 0, "data",      data     );
  if(has_variance)
    mxSetField(plhs[0], 0, "variance",  variance );
  mxSetField(plhs[0], 0, "time",      time     );
  mxSetField(plhs[0], 0, "label",     label    );

  free(dat);
  if(has_variance)
    free(var);
 
  /* close the data file */
  eep_free(hdr);
  eepio_fclose(fp);

  return;
}
