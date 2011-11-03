#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2009-2011 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4541 2011-11-03 16:37:35Z guillaume $

include ../src/Makefile.var

all: fieldtrip-all bemcp-all peer-all

clean: fieldtrip-clean bemcp-clean peer-clean

distclean: clean fieldtrip-distclean bemcp-distclean peer-distclean

install: fieldtrip-install bemcp-install peer-instal

tarball: fieldtrip-tarball bemcp-tarball peer-tarball
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
        
###############################################################################
# PEER
###############################################################################
.PHONY: peer
peer peer-%:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        In peer"
	@ echo "_____________________________________________________________"
	@ echo ""
	@ $(MAKE) -f Makefile.peer --no-print-directory $*;
