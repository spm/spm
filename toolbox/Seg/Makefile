#!/usr/bin/env make -f
# Seg Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4037 2010-08-06 14:52:27Z guillaume $

include ../../src/Makefile.var

SPMMEX = spm_mrf.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
