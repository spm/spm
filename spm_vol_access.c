#ifndef lint
static char sccsid[] = "%W% John Ashburner %E%";
#endif

#define UNSIGNED_CHAR   2
#define SIGNED_SHORT    4
#define SIGNED_INT      8
#define FLOAT   16
#define DOUBLE  64

#define SIGNED_SHORT_S	SIGNED_SHORT<<8
#define SIGNED_INT_S	SIGNED_INT<<8
#define FLOAT_S	FLOAT<<8
#define DOUBLE_S	DOUBLE<<8

void resample(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset, dtype)
int m, xdim,ydim,zdim, hold, dtype;
double out[], x[], y[], z[], background, scale,offset;
unsigned char vol[];
{
	extern void resample_uchar(), resample_short(), resample_int(), resample_float(),
		resample_double(), resample_short_s(), resample_int_s(), resample_float_s(),
		resample_double_s();

	if (dtype == UNSIGNED_CHAR)
		 resample_uchar(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_SHORT)
		 resample_short(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_INT)
		   resample_int(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == FLOAT)
		 resample_float(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == DOUBLE)
		resample_double(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_SHORT_S)
		resample_short_s(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_INT_S)
		resample_int_s(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == FLOAT_S)
		resample_float_s(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == DOUBLE_S)
		resample_double_s(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else
		exit(1);
}

void resample_d(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset, dtype)
int m, xdim,ydim,zdim, hold, dtype;
double out[], gradx[],grady[],gradz[], x[], y[], z[], background, scale,offset;
unsigned char vol[];
{
	extern void resample_d_uchar(), resample_d_short(), resample_d_int(), resample_d_float(),
		resample_d_double(), resample_d_short_s(), resample_d_int_s(), resample_d_float_s(),
		resample_d_double_s();

	if (dtype == UNSIGNED_CHAR)
		 resample_d_uchar(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_SHORT)
		 resample_d_short(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_INT)
		   resample_d_int(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == FLOAT)
		 resample_d_float(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == DOUBLE)
		resample_d_double(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_SHORT_S)
		resample_d_short_s(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == SIGNED_INT_S)
		resample_d_int_s(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == FLOAT_S)
		resample_d_float_s(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else if (dtype == DOUBLE_S)
		resample_d_double_s(m,vol,out,gradx,grady,gradz,x,y,z,xdim,ydim,zdim, hold, background, scale,offset);
	else
		exit(1);
}

slice(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset, dtype)
int ydim1,xdim1, xdim2,ydim2,zdim2, hold, dtype;
double image[], mat[], background, scale,offset;
unsigned char vol[];
{
	int sts = 1;
	if (dtype == UNSIGNED_CHAR)
		 sts = slice_uchar(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == SIGNED_SHORT)
		 sts = slice_short(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == SIGNED_INT)
		 sts = slice_int(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == FLOAT)
		sts = slice_float(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == DOUBLE)
		sts = slice_double(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == SIGNED_SHORT_S)
		sts = slice_short_s(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == SIGNED_INT_S)
		sts = slice_int_s(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == FLOAT_S)
		sts = slice_float_s(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else if (dtype == DOUBLE_S)
		sts = slice_double_s(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale,offset);
	else
		sts = 1;
	return(sts);
}

