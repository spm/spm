#!/usr/bin/env make -f
# FieldMap Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 6791 2016-04-28 14:47:20Z john $

include ../../src/Makefile.var

SPMMEX =\
	pm_invert_phasemap_dtj.$(SUF) pm_merge_regions.$(SUF) \
	pm_create_connectogram_dtj.$(SUF) pm_pad.$(SUF) \
	pm_estimate_ramp.$(SUF) pm_restore_ramp.$(SUF) \
	pm_ff_unwrap.$(SUF) 

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
