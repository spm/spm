#!/usr/bin/env make -f
# XMLTree Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 7532 2019-02-14 12:03:24Z guillaume $

include ../../src/Makefile.var

SPMMEX = xml_findstr.$(MEXEXT)

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
