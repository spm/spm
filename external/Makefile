#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 5914 2014-03-12 12:17:21Z guillaume $

include ../src/Makefile.var

all: fieldtrip-all bemcp-all

clean: fieldtrip-clean bemcp-clean

distclean: clean fieldtrip-distclean bemcp-distclean

install: fieldtrip-install bemcp-install

tarball: fieldtrip-tarball bemcp-tarball
	$(TAR) Af spm_mex.tar fieldtrip_mex.tar
	$(TAR) Af spm_mex.tar bemcp/spm_mex.tar

###############################################################################
# FIELDTRIP
###############################################################################
.PHONY: fieldtrip
fieldtrip fieldtrip-all fieldtrip-clean fieldtrip-distclean fieldtrip-install fieldtrip-tarball:
	@ $(call verb,fieldtrip)
	@ $(MAKE) -f Makefile.fieldtrip --no-print-directory $(subst fieldtrip,,$(subst -,,$@))

###############################################################################
# BEMCP
###############################################################################
.PHONY: bemcp
bemcp bemcp-all bemcp-clean bemcp-distclean bemcp-install bemcp-tarball:
	@ $(call verb,bemcp)
	@ $(MAKE) -C bemcp --no-print-directory $(subst bemcp,,$(subst -,,$@))
