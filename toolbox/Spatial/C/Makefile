include Makefile.var

#PTXFILES = pushpull.ptx lbessi.ptx operator.ptx shootfun.ptx TVdenoise3d.ptx TVdenoise2d.ptx TVdenoise2d_old.ptx
PTXFILES = pushpull.ptx lbessi.ptx operator.ptx TVdenoise3d.ptx TVdenoise2d.ptx
LIBFILES = pushpull.$(SOSUF) TVdenoise3d.$(SOSUF)

#CC       = x86_64-w64-mingw32-gcc-win32

all: $(PTXFILES) $(LIBFILES)

install: all $(PTXDIR) $(LIBDIR)
	$(COPY) $(PTXFILES) $(PTXDIR)
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

%.so : %.c
	$(CC) -O4 -shared -fPIC -Wl,-soname,$@ -o $@ -lm $<

%.dll : %.c
	$(CC) -O4 -shared -fPIC -Wl,-soname,$@ -o $@ -lm $<

#%.so : %.c
#	gcc -O4 -shared -fPIC -Wl,-soname,$@ -o $@ -lm $<

%.ptx : %.cu
	$(NVCC) -ptx --ptxas-options --verbose $< 

$(PTXDIR) :
	$(MKDIR) $(PTXDIR)

$(LIBDIR) :
	$(MKDIR) $(LIBDIR)



