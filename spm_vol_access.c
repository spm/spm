#ifndef lint
static char sccsid[] = "%W% John Ashburner %E%";
#endif
#include "spm_vol_utils.h"


void resample(m,vol,out,x,y,z,hold, background)
int m, hold;
double out[], x[], y[], z[], background;
MAPTYPE *vol;
{
	extern void resample_uchar(), resample_short(), resample_int(), resample_float(),
		resample_double(), resample_short_s(), resample_int_s(), resample_float_s(),
		resample_double_s();

	if (vol->dtype == UNSIGNED_CHAR)
		 resample_uchar(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_SHORT)
		 resample_short(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_INT)
		   resample_int(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == FLOAT)
		 resample_float(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == DOUBLE)
		resample_double(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_SHORT_S)
		resample_short_s(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_INT_S)
		resample_int_s(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == FLOAT_S)
		resample_float_s(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == DOUBLE_S)
		resample_double_s(m,vol->data,out,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else
		exit(1);
}

void resample_d(m,vol,out,gradx,grady,gradz,x,y,z, hold, background)
int m, hold;
double out[], gradx[],grady[],gradz[], x[], y[], z[], background;
MAPTYPE *vol;
{
	extern void resample_d_uchar(), resample_d_short(), resample_d_int(), resample_d_float(),
		resample_d_double(), resample_d_short_s(), resample_d_int_s(), resample_d_float_s(),
		resample_d_double_s();

	if (vol->dtype == UNSIGNED_CHAR)
		 resample_d_uchar(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_SHORT)
		 resample_d_short(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_INT)
		   resample_d_int(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == FLOAT)
		 resample_d_float(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == DOUBLE)
		resample_d_double(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_SHORT_S)
		resample_d_short_s(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_INT_S)
		resample_d_int_s(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == FLOAT_S)
		resample_d_float_s(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else if (vol->dtype == DOUBLE_S)
		resample_d_double_s(m,vol->data,out,gradx,grady,gradz,x,y,z,vol->dim[0],vol->dim[1],vol->dim[2],
			hold, background, vol->scale,vol->offset);
	else
		exit(1);
}

slice(mat, image, xdim1,ydim1, vol, hold,background)
int ydim1,xdim1, hold;
double image[], mat[], background;
MAPTYPE *vol;
{
	int sts = 1;
	if (vol->dtype == UNSIGNED_CHAR)
		 sts = slice_uchar(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2],
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_SHORT)
		 sts = slice_short(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2],
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_INT)
		 sts = slice_int(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2],
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == FLOAT)
		sts = slice_float(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2],
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == DOUBLE)
		sts = slice_double(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2],
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_SHORT_S)
		sts = slice_short_s(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2], 
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == SIGNED_INT_S)
		sts = slice_int_s(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2], 
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == FLOAT_S)
		sts = slice_float_s(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2],
			hold,background, vol->scale,vol->offset);
	else if (vol->dtype == DOUBLE_S)
		sts = slice_double_s(mat, image, xdim1,ydim1, vol->data, vol->dim[0],vol->dim[1],vol->dim[2], 
			hold,background, vol->scale,vol->offset);
	else
		sts = 1;
	return(sts);
}

