/* spm_map.h */

#include <sys/types.h>
#ifdef SPM_WIN32
#include "win32mmap.h"
#endif

#define MAGIC 110494

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

typedef struct map
{
	double dim[3];
	double pixdim[3];
	double scale;
	int dtype;
	int off;
	caddr_t map;
	unsigned char *data;
	size_t len;
	int magic;
	int prot;
	int flags;
	pid_t pid;
}	MAP;
