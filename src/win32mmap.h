/* %W% %E%
mmap replacement for 32 bit windows
does not implement mmap, but map_file 
See win32mmap.c for comments
*/

#ifndef WIN32MMAP_DEF
#define WIN32MMAP_DEF

#include <sys/types.h>

/* mmap constants */
#define PROT_READ	1
#define PROT_WRITE	2
#define MAP_SHARED	1
#define MAP_PRIVATE 2
#define MAP_FIXED	4
#define MAP_NORESERVE	8

/* not defined in mingw headers */
typedef char* caddr_t;

/* function prototypes - see map_file.c for comments*/
caddr_t map_file(char *filename, caddr_t addr, size_t length,
		int prot, int flags, off_t offset);

int unmap_file(caddr_t start);

#endif /* WIN32MMAP_DEF */

