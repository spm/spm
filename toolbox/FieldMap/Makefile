#!/usr/bin/env make -f
# FieldMap Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 2457 2008-11-12 09:24:30Z guillaume $

include ../../src/Makefile.var

SPMMEX =\
	pm_invert_phasemap_dtj.$(SUF) pm_merge_regions.$(SUF) \
	pm_create_connectogram_dtj.$(SUF) pm_pad.$(SUF) \
	pm_estimate_ramp.$(SUF) pm_restore_ramp.$(SUF) \
	pm_ff_unwrap.$(SUF) pm_smooth_phasemap_dtj.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
