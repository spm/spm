#!/usr/bin/env make -f
# HDW Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 1288 2008-04-01 16:05:52Z guillaume $

include ../../src/Makefile.var

SPMMEX =\
	spm_affdef.$(SUF) spm_def2det.$(SUF) spm_warp.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
	@ $(CHMODIT) $@
