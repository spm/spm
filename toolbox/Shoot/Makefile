#!/usr/bin/env make -f
# Shoot Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4588 2011-12-08 18:28:59Z guillaume $

include ../../src/Makefile.var

SPMMEX = shoot3.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

shoot3.$(SUF): shoot3.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c shoot_regularisers.c shoot_expm3.c shoot_dartel.c shoot_boundary.c shoot_boundary.h  shoot_dartel.h  shoot_diffeo3d.h  shoot_expm3.h  shoot_multiscale.h  shoot_optim3d.h  shoot_regularisers.h
	$(MEX) shoot3.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c shoot_regularisers.c shoot_expm3.c shoot_dartel.c shoot_boundary.c $(MEXEND)
