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

%.$(SOSUF) : %.c
	$(CC) -O3 -shared -fPIC -Wl,-$(SONAME),$@ -o $@ -lm $<

%.ptx : %.cu
	$(NVCC) -ptx --ptxas-options --verbose $< 

$(PTXDIR) :
	$(MKDIR) $(PTXDIR)

$(LIBDIR) :
	$(MKDIR) $(LIBDIR)
