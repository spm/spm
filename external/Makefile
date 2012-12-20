#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 5148 2012-12-20 18:42:30Z guillaume $

include ../src/Makefile.var

all: fieldtrip-all bemcp-all #engine-all

clean: fieldtrip-clean bemcp-clean #engine-clean

distclean: clean fieldtrip-distclean bemcp-distclean #engine-distclean

install: fieldtrip-install bemcp-install #engine-install

tarball: fieldtrip-tarball bemcp-tarball #engine-tarball
	$(TAR) Af spm_mex.tar fieldtrip_mex.tar
	$(TAR) Af spm_mex.tar bemcp/spm_mex.tar
#	$(TAR) Af spm_mex.tar engine_mex.tar

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
        
###############################################################################
# ENGINE
###############################################################################
.PHONY: engine
engine engine-%:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        In engine"
	@ echo "_____________________________________________________________"
	@ echo ""
	@ $(MAKE) -f Makefile.engine --no-print-directory $*;
