#!/usr/bin/env make -f
# EEG-BEM Makefile called by {SPM}/external/Makefile
#
# Christophe Phillips
# $Id: Makefile 7532 2019-02-14 12:03:24Z guillaume $

include ../../src/Makefile.var

BEMCPMEX = \
	bem_Cii_cog.$(MEXEXT) bem_Cii_cst.$(MEXEXT) bem_Cii_lin.$(MEXEXT) \
	bem_Cij_cog.$(MEXEXT) bem_Cij_cst.$(MEXEXT) bem_Cij_lin.$(MEXEXT) \
	bem_Gi_cog.$(MEXEXT) bem_Gi_vert.$(MEXEXT)

all: $(BEMCPMEX)

clean:

distclean: clean
	$(DEL) $(BEMCPMEX)

install:

tarball: all
	$(TAR) -cf spm_mex.tar --transform='s,^,bemcp/,g' $(BEMCPMEX)

%.$(MEXEXT) : %.c
	$(MEX) $< $(MEXEND)
