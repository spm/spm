/* win32mmap.c 
Windows 32bit partial emulation of Unix mmap, munmap functions

mmap -> map_file
munmap -> unmap_file
 
map_file differs from mmap in that
1) it accepts a filename, not a file descriptor, because the Win32
	memory mapping requires a windows file handle
2) it does half an attempt to implement the flags parameter of mmap, 
	because there are not straightforward Win32 equivalents

  mmap has the following syntax for comparison:
  caddr_t mmap(caddr_t addr, size_t len, int prot, int flags,
          int fildes, off_t off);

unmap_file is very similar to munmap, except it does not
want the (second) length parameter, which is not easily
implemented using Win32 
*/
#ifndef lint
static char sccsid[] = "@(#)win32mmap.c	2.1 Matthew Brett 99/01/15";
#endif

#include <windows.h>
#include "win32mmap.h"

caddr_t map_file(char *filename, caddr_t addr, size_t length,
		int prot, int flags, off_t offset) {
  HANDLE file_handle, file_mapping;
  caddr_t retval, staddr;
  int iAccess, iProt, iFlags;

  /* default return is addr -1 (as for mmap) */
  retval = addr - 1;
  
  /* recode PROT and FLAGS into desired access */
  iAccess = 0;
  iProt = 0;
  iFlags = 0;
  if (prot & PROT_READ) { 
	iAccess = iAccess | GENERIC_READ;
	iProt = PAGE_READONLY;		/* READONLY assumed */
	iFlags = FILE_MAP_READ;	/* and FLAGS are ignored */
	}
  if (prot & PROT_WRITE) {
	iAccess = iAccess | GENERIC_WRITE;
	if (flags & MAP_PRIVATE)  {	/* Write changes won't be shared */
		iProt = PAGE_WRITECOPY;
		iFlags = FILE_MAP_COPY;
	} 
	else if (flags & MAP_SHARED) { /* Write changes will be shared */
		iProt = PAGE_READWRITE;
		iFlags = FILE_MAP_WRITE;
	}
  }
  if ((iProt != 0) & (iFlags != 0)) {	/* valid PROT, FLAGS options passed */
		
  	  file_handle = CreateFile(filename,
				   iAccess,
				   FILE_SHARE_WRITE,
				   NULL,
				   OPEN_EXISTING,
				   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS,
				   NULL);
	  if (file_handle == NULL) 
		  return retval;	/* "Could not open file"*/
	  file_mapping = CreateFileMapping(file_handle,
					   NULL,
					   iProt,
					   0, 0, 
					   NULL);

	  if (file_mapping == NULL)
		  return retval;	/* "Could not create file mapping" */
	  staddr = (caddr_t)MapViewOfFileEx(file_mapping,
				  iFlags,
				  0, offset, length, (LPVOID) addr);
	  if (staddr == NULL)
		  return retval;	/*"Could not map the file"*/
	  if (!CloseHandle(file_mapping))
		  return retval;	/*"Could not close file mapping"*/
	  if (!CloseHandle(file_handle))
		  return retval;	/*"Could not close file handle"*/
	  retval = staddr;		/* return address of mapped memory */
	}
	return (caddr_t) retval;
}

int unmap_file(caddr_t start){
	return (int)(UnmapViewOfFile((LPVOID) start));
}
