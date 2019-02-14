#!/usr/bin/env make -f
# GIfTI Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 7532 2019-02-14 12:03:24Z guillaume $

include ../../src/Makefile.var

SPMMEX = zstream.$(MEXEXT)

ifeq (mex,$(MEXEXT))
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

%.$(MEXEXT) : %.c
	$(MEX) $< $(MEXEND)
