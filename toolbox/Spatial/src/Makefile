# The lib files can be compiled through github actions
# but this is not done for the ptx files as it requires
# payment. Therefore some messing about is needed in terms
# of downloading the artifacts from github (with potential
# problems due to security), unzipping and making ptx
# ptx fies locally. The tar.gz file is then generated and
# the Artifacts.toml edited to contain the appropriate SHA.

include Makefile.var

PTXFILES = pushpull.ptx lbessi.ptx sparse_operator.ptx TVdenoise3d.ptx TVdenoise2d.ptx blip.ptx
LIBFILES = pushpull.$(SOSUF) TVdenoise3d.$(SOSUF) sparse_operator.$(SOSUF) blip.$(SOSUF)
PPLIB    = ../pp_lib.tar.gz

all: $(LIBFILES)

ptx: $(PTXFILES)

install: all $(LIBDIR)
	$(COPY) $(LIBFILES) $(LIBDIR)

install-ptx: ptx $(PTXDIR) 
	$(COPY) $(PTXFILES) $(PTXDIR)

tarfile:
	tar -C .. -cvf - lib | gzip -v > $(PPLIB)
	julia get_sha.jl $(PPLIB)

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



