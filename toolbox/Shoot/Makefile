#!/usr/bin/env make -f
# Shoot Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4871 2012-08-30 14:11:53Z john $

include ../../src/Makefile.var

SPMMEX = shoot3.$(SUF) optimN_mex.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

shoot3.$(SUF): shoot3.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c shoot_regularisers.c shoot_expm3.c shoot_invdef.c shoot_dartel.c shoot_boundary.c shoot_bsplines.c ../../src/bsplines.c
	$(MEX) shoot3.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c shoot_regularisers.c shoot_expm3.c shoot_invdef.c shoot_dartel.c shoot_boundary.c shoot_bsplines.c ../../src/bsplines.c -I../../src/ -DIMAGE_SINGLE $(MEXEND)

optimN_mex.$(SUF): optimN_mex.c  shoot_optimN.c shoot_multiscale.c shoot_boundary.c
	$(MEX)  optimN_mex.c shoot_optimN.c shoot_multiscale.c shoot_boundary.c $(MEXEND)
