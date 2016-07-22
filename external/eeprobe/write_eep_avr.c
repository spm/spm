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
#include <eep/eepmisc.h>		/* MPI-ANT specific */
#include <cnt/cnt.h>		/* MPI-ANT specific */

#define NAME "write_eep_avr"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  /* these variables are MPI-ANT specific for reading the data */
  char filename[256];
  char lab[8];
  int nchan;
  int nsample;
  int trialc;
  int nsweeps;
  char condlab[256];
  char condcol[256];
  float psi=0;
  float rate;
  double *mx_data_ptr,
         *mx_variance_ptr;
  FILE *fp;
  float *buf;
  eeg_t *avr;
  eegchan_t *channel_structure;
  int i,j;

  /* these variables are Matlab specific for interfacing */
  mxArray *mx_data;
  mxArray *mx_variance;
  mxArray *mx_label;
  mxArray *mx_tmp;

  /*****************************
   * test proper function call *
   *****************************/
  if (nrhs!=2)
    mexErrMsgTxt ("Invalid number of input arguments");

  /**************************************************************************************
   * be sure to setup important variables here.
   * stuff like number of channels(nchan), the names for the channels(in label), etc... *
   **************************************************************************************/
  nchan = (int)*(mxGetPr(mxGetField(prhs[1], 0, "nchan")));
  nsample = (int)*(mxGetPr(mxGetField(prhs[1], 0, "npnt")));
  rate = (float)(*mxGetPr(mxGetField(prhs[1], 0, "rate")));
  trialc = (int)*(mxGetPr(mxGetField(prhs[1], 0, "trialc")));
  nsweeps = (int)*(mxGetPr(mxGetField(prhs[1], 0, "nsweeps")));
  mxGetString(mxGetField(prhs[1], 0, "condlab"), condlab, 256);
  mxGetString(mxGetField(prhs[1], 0, "condcol"), condcol, 256);
  if(mxGetField(prhs[1], 0, "psi")!=NULL) {
    psi = (float)(*mxGetPr(mxGetField(prhs[1], 0, "psi")));
  }
  mx_label = mxGetField(prhs[1], 0, "label");
  mx_data=mxGetField(prhs[1], 0, "data");
  mx_variance=mxGetField(prhs[1], 0, "variance");
  mx_data_ptr=mxGetPr(mx_data);
  if(mx_variance!=NULL) {
    mx_variance_ptr=mxGetPr(mx_variance);
  }

  /******************
   * build channels *
   ******************/
  channel_structure=eep_chan_init(nchan);
  if(channel_structure==NULL) {
    mexErrMsgTxt ("could not init channels");
  }
  for(i=0;i<nchan;i++) {
    mxGetString(mxGetCell(mx_label, i), lab, 8);
    eep_chan_set(channel_structure, i, lab , 1, 0.001, "uV");
  }

  /******************
   * initialize avr *
   ******************/
  avr=eep_init_from_values(1.0f/rate, nchan, channel_structure);
  if(channel_structure==NULL) {
    mexErrMsgTxt ("could not init avr file");
  }

  eep_set_total_trials(avr, trialc);
  eep_set_averaged_trials(avr, nsweeps);
  eep_set_conditionlabel(avr, condlab);
  eep_set_conditioncolor(avr, condcol);
  eep_set_pre_stimulus_interval(avr, psi);

  /*************
   * open file *
   *************/
  mxGetString(prhs[0], filename, 256);
  fp=eepio_fopen(filename,"wb");
  if(fp==NULL) {
    mexErrMsgTxt ("Could not open file");
  }
  if(eep_create_file(avr, filename, fp, NULL, 0, "avr matlab exporter")) {
    mexErrMsgTxt ("Could not create avr file");
  }
  /**************************************
   * setup libeep to start writing data *
   **************************************/
  if(eep_prepare_to_write(avr, DATATYPE_AVERAGE, rate, NULL) != CNTERR_NONE) {
    mexErrMsgTxt ("Could not setup avr for writing data");
  }
  /*******************
   * allocate buffer *
   *******************/
  buf=(float*)(malloc(CNTBUF_SIZE(avr, 1)));
  if(buf==NULL) {
    mexErrMsgTxt ("Could not alocate buffer");
  }
  /****************************
   * iterate through data ... *
   ****************************/
  for(i=0;i<nsample;i++) {
    /* build buffer */
    for(j=0;j<nchan;j++) {
      buf[j] = mx_data_ptr[i*nchan+j] * 1000.0;
    }
    /* write */
    if(eep_write_float(avr, buf, 1) != CNTERR_NONE) {
      mexErrMsgTxt ("Could not write sample");
    }
  }
  /********************
   * ... and variance *
   ********************/
  if(mx_variance!=NULL) {
    if(eep_prepare_to_write(avr, DATATYPE_STDDEV, rate, NULL) != CNTERR_NONE) {
      mexErrMsgTxt ("Could not prepar variance data");
    }
    for(i=0;i<nsample;i++) {
      /* build buffer */
      for(j=0;j<nchan;j++) {
        buf[j] = mx_variance_ptr[i*nchan+j] * 1000.0;
      }
      /* write */
      if(eep_write_float(avr, buf, 1) != CNTERR_NONE) {
        mexErrMsgTxt ("Could not write sample");
      }
    }
  }
  /*************************
   * Done! shutdown nicely *
   *************************/
  eep_finish_file(avr);
  free(buf);
  eepio_fclose(fp);

  return;
}
