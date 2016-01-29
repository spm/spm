#!/usr/bin/env make -f
# GIfTI Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 6704 2016-01-29 17:23:13Z guillaume $

include ../../src/Makefile.var

SPMMEX = zstream.$(SUF)

ifeq (mex,$(SUF))
  export CFLAGS = $(shell $(MEX) -p CFLAGS) -std=c99
else
  ifeq (windows,$(PLATFORM))
    MEXOPTS += CFLAGS=\"$$CFLAGS -std=c99\"
  else
    MEXOPTS += CFLAGS='$$CFLAGS -std=c99'
  endif
endif

all: $(SPMMEX)
	@:

clean:
	@:
        
distclean: clean
	$(DEL) $(SPMMEX)

install:
	@:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
