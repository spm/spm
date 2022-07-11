#!/usr/bin/env make -f
# @file_array Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

include ../../src/Makefile.var

SPMMEX = file2mat.$(MEXEXT) mat2file.$(MEXEXT) init.$(MEXEXT)

all: $(SPMMEX)
	@:

clean:
	@:
ifeq (mex,$(MEXEXT))
	$(DEL) $(subst .$(MEXEXT),.o,$(SPMMEX))
endif

distclean: clean
	$(DEL) $(SPMMEX)

install:
	@:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(MEXEXT) : %.c
	$(MEX) $< $(MEXEND)
