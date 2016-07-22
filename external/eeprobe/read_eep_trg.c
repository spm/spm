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
#include <cnt/trg.h>		/* MPI-ANT specific */

#define NAME "read_eep_trg"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  /* these variables are MPI-ANT specific for reading the data */
  char filename[256];
  int status;
  FILE *trigger_file;
  trg_t* trigger_info;
  double period;
  short chanc;
  int i;

  /* these variables are Matlab specific for interfacing */
  mxArray *data;

  const int nfields = 4;
  const char *field_names[] = {
    "time",
    "offset",
    "code",
    "type"};

  if (nrhs!=1)
    mexErrMsgTxt ("Invalid number of input arguments");

  mxGetString(prhs[0], filename, 256);

  /* open the data file */
 if ((trigger_file = eepio_fopen(filename, "r"))==NULL)
    mexErrMsgTxt ("Could not open file");

  trigger_info = trg_file_read_unchecked(trigger_file, &period, &chanc);

  /* create the struct array with dimensions 1x1 */
  data = mxCreateStructMatrix(1, trigger_info->c, nfields, field_names);


  for(i=0; i< trigger_info->c; i++)
  {
    trgentry_t* entry = &(trigger_info->v[i]);
    mxArray * time;
    mxArray * offset;
    mxArray * code;
    mxArray * type;

    time = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(time) = (double) entry->sample * period;

    offset = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(offset) = (double) entry->sample+1;

    code = mxCreateString(entry->code);

    type = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(type) = (double) atof(entry->code);

    /* fill the struct array with the variables */
    mxSetField(data, i, "time", time );
    mxSetField(data, i, "offset", offset );
    mxSetField(data, i, "code", code );
    mxSetField(data, i, "type", type );

  }
  eepio_fclose(trigger_file);

  plhs[0]  = data;

  return;
}
