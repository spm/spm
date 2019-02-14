#!/usr/bin/env make -f
# FieldMap Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 7532 2019-02-14 12:03:24Z guillaume $

include ../../src/Makefile.var

SPMMEX =\
	pm_invert_phasemap_dtj.$(MEXEXT) pm_merge_regions.$(MEXEXT) \
	pm_create_connectogram_dtj.$(MEXEXT) pm_pad.$(MEXEXT) \
	pm_estimate_ramp.$(MEXEXT) pm_restore_ramp.$(MEXEXT) \
	pm_ff_unwrap.$(MEXEXT) 

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
