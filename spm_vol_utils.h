#define UNSIGNED_CHAR     2
#define SIGNED_SHORT      4
#define SIGNED_INT        8
#define FLOAT             16
#define DOUBLE            64
#define SIGNED_SHORT_S    (SIGNED_SHORT<<8)
#define SIGNED_INT_S      (SIGNED_INT<<8)
#define FLOAT_S           (FLOAT<<8)
#define DOUBLE_S          (DOUBLE<<8)
#define SIGNED_CHAR       (UNSIGNED_CHAR+128) 
#define UNSIGNED_SHORT    (SIGNED_SHORT+128) 
#define UNSIGNED_INT      (SIGNED_INT+128)
#define UNSIGNED_SHORT_S  (UNSIGNED_SHORT<<8)
#define UNSIGNED_INT_S    (UNSIGNED_INT<<8)

#include <sys/types.h>
#include <sys/mman.h>

typedef struct maptype
{
	int dim[3];		/* Dimensions of the volume */
	double *scale, *offset;	/* Scalefactor and offset, such that true_intensity = vox*scale+offset */
	int dtype;		/* Data-type of volume */
	void **data;	/* Pointer to data */
	double mat[16];

	caddr_t addr;
	size_t len;
}	MAPTYPE;
