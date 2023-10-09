include Makefile.var

#PTXFILES = pushpull.ptx lbessi.ptx operator.ptx shootfun.ptx TVdenoise3d.ptx TVdenoise2d.ptx TVdenoise2d_old.ptx
PTXFILES = pushpull.ptx lbessi.ptx operator.ptx TVdenoise3d.ptx TVdenoise2d.ptx
LIBFILES = pushpull.$(SOSUF) TVdenoise3d.$(SOSUF)

CUDA_VERSION := $(shell $(NVCC) --version 2> /dev/null)
ifndef CUDA_VERSION
	PTXFILES =
	PTXDIR =
endif

all: $(PTXFILES) $(LIBFILES)

install: all $(PTXDIR) $(LIBDIR)
ifdef CUDA_VERSION
	$(COPY) $(PTXFILES) $(PTXDIR)
endif
	$(COPY) $(LIBFILES) $(LIBDIR)

clean:
	$(DEL) $(PTXFILES) $(LIBFILES)

#pushpull.$EXT : pushpull.c pushpull_dev.cu
#	$(CC) -Wall -fPIC -O3 -c pushpull.c
#	$(CC) -shared -Wl,-soname,pushpull.$EXT -o pushpull.$EXT pushpull.o

#TVdenoise3d.$EXT : TVdenoise3d.c TVdenoise3d_dev.cu
#	$(CC) -Wall -fPIC -O3 -c TVdenoise3d.c
#	$(CC) -shared -Wl,-soname,TVdenoise3d.$EXT -o TVdenoise3d.$EXT TVdenoise3d.o

#pushpull.so : pushpull.c shootfun.c
#	cc -mwin32 -O4 -shared -fPIC -Wl,-soname,pushpull.so -o pushpull.so -lm pushpull.c shootfun.c

#pushpull.so : pushpull.c
#	$(CC) -Wall -fPIC -O3 -c pushpull.c
#	$(CC) -shared -Wl,-soname,pushpull.so -o pushpull.so pushpull.o

#TVdenoise3.so : TVdenoise3.c
#	$(CC) -Wall -fPIC -O3 -c TVdenoise3.c
#	$(CC) -shared -Wl,-soname,TVdenoise3.so -o TVdenoise3.so TVdenoise3.o

%.$(SOSUF) : %.c
	$(CC) -O4 -shared -fPIC -Wl,-$(SONAME),$@ -o $@ -lm $<

%.ptx : %.cu
	$(NVCC) -ptx --ptxas-options --verbose $< 

$(PTXDIR) :
	$(MKDIR) $(PTXDIR)

$(LIBDIR) :
	$(MKDIR) $(LIBDIR)
