#!/usr/bin/env make -f
# @file_array Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 5446 2013-04-24 16:56:51Z guillaume $

include ../../src/Makefile.var

SPMMEX = file2mat.$(SUF) mat2file.$(SUF) init.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
