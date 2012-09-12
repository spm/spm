#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4918 2012-09-12 15:28:50Z guillaume $

include ../src/Makefile.var

all: fieldtrip-all bemcp-all

clean: fieldtrip-clean bemcp-clean

distclean: clean fieldtrip-distclean bemcp-distclean

install: fieldtrip-install bemcp-install

tarball: fieldtrip-tarball bemcp-tarball
	$(TAR) Af spm_mex.tar fieldtrip_mex.tar
	$(TAR) Af spm_mex.tar bemcp/spm_mex.tar
	$(TAR) Af spm_mex.tar peer_mex.tar

###############################################################################
# FIELDTRIP
###############################################################################
.PHONY: fieldtrip
fieldtrip fieldtrip-%:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        In fieldtrip"
	@ echo "_____________________________________________________________"
	@ echo ""
	@ $(MAKE) -f Makefile.fieldtrip --no-print-directory $*;

###############################################################################
# BEMCP
###############################################################################
.PHONY: bemcp
bemcp bemcp-%:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        In bemcp"
	@ echo "_____________________________________________________________"
	@ echo ""
	@ $(MAKE) -C bemcp --no-print-directory $*;
