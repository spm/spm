#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 5182 2013-01-09 12:00:59Z guillaume $

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
        
###############################################################################
# ENGINE
###############################################################################
.PHONY: engine
engine engine-all engine-clean engine-distclean engine-install engine-tarball:
	@ $(call verb,engine)
	@ $(MAKE) -f Makefile.engine --no-print-directory $(subst engine,,$(subst -,,$@))

###############################################################################
# MESSAGE
###############################################################################
define verb
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        In " $1
	@ echo "_____________________________________________________________"
	@ echo ""
endef
