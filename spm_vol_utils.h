#define UNSIGNED_CHAR   2
#define SIGNED_SHORT    4
#define SIGNED_INT      8
#define FLOAT           16
#define DOUBLE          64
#define SIGNED_SHORT_S  SIGNED_SHORT<<8
#define SIGNED_INT_S    SIGNED_INT<<8
#define FLOAT_S         FLOAT<<8
#define DOUBLE_S        DOUBLE<<8

typedef struct maptype
{
	int dim[3];		/* Dimensions of the volume */
	double pixdim[3];	/* Voxel sizes of the volume */
	double scale, offset;	/* Scalefactor and offset, such that true_intensity = vox*scale+offset */
	int dtype;		/* Data-type of volume */
	unsigned char *data;	/* Pointer to data */
}	MAPTYPE;
