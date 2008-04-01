#!/usr/bin/env make -f
# DARTEL Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 1287 2008-04-01 14:50:25Z guillaume $

SPMMEX =\
	dartel2.$(SUF) dartel3.$(SUF) \
	optim1.$(SUF) optimN.$(SUF) optimNn.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

dartel2.$(SUF): dartel2.c optimizer2d.c diffeo2d.c
	$(MEX) dartel2.c optimizer2d.c diffeo2d.c $(MEXEND)
	@ $(CHMODIT) $@

dartel3.$(SUF): dartel3.c optimizer3d.c diffeo3d.c
	$(MEX) dartel3.c optimizer3d.c diffeo3d.c $(MEXEND)
	@ $(CHMODIT) $@

optim1.$(SUF): optim1.c optim1_mex.c 
	$(MEX) optim1.c optim1_mex.c $(MEXEND)
	@ $(CHMODIT) $@

optimN.$(SUF): optimN.c optimN_mex.c 
	$(MEX) optimN.c optimN_mex.c $(MEXEND)
	@ $(CHMODIT) $@

optimNn.$(SUF): optimN.c optimN_mex.c
	$(MEX) optimN.c optimN_mex.c -DNEUMANN -output optimNn $(MEXEND)
	@ $(CHMODIT) $@
