#!/usr/bin/env make -f
# EEG-BEM Makefile called by {SPM}/src/Makefile
#
# Christophe Phillips
# $Id: Makefile 2718 2009-02-09 19:14:20Z vladimir $

include ../../src/Makefile.var

SPMMEX = \
	bem_Cii_cog.$(SUF) bem_Cii_cst.$(SUF) bem_Cii_lin.$(SUF) \
	bem_Cij_cog.$(SUF) bem_Cij_cst.$(SUF) bem_Cij_lin.$(SUF) \
	bem_Gi_cog.$(SUF) bem_Gi_vert.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
