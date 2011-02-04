#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4193 2011-02-04 13:51:29Z guillaume $

include ../src/Makefile.var

FIELDTRIPMEX =\
	fieldtrip/src/read_24bit.$(SUF) \
	fieldtrip/src/mxSerialize.$(SUF) \
	fieldtrip/src/mxDeserialize.$(SUF) fieldtrip/src/read_16bit.$(SUF) \
	fieldtrip/fileio/@uint64/max.$(SUF) fieldtrip/fileio/@uint64/min.$(SUF) \
	fieldtrip/fileio/@uint64/minus.$(SUF) fieldtrip/fileio/@uint64/plus.$(SUF) \
	fieldtrip/fileio/@uint64/rdivide.$(SUF) fieldtrip/fileio/@uint64/times.$(SUF) \
	fieldtrip/fileio/@uint64/abs.$(SUF) \
	\
	fieldtrip/src/meg_leadfield1.$(SUF) fieldtrip/src/ptriproj.$(SUF) \
	fieldtrip/src/lmoutr.$(SUF) fieldtrip/src/plgndr.$(SUF) \
	fieldtrip/src/solid_angle.$(SUF) fieldtrip/src/routlm.$(SUF) \
	\
	fieldtrip/src/solid_angle.$(SUF) \
	\
	fieldtrip/src/ltrisect.$(SUF) fieldtrip/src/plinproj.$(SUF) \
	fieldtrip/src/splint_gh.$(SUF) \
	\
	fieldtrip/src/read_ctf_shm.$(SUF) fieldtrip/src/write_ctf_shm.$(SUF) # Unix only
	
BEMCP = bemcp

all: $(FIELDTRIPMEX) bemcp-all

clean: bemcp-clean

distclean: clean bemcp-distclean
	$(DEL) $(FIELDTRIPMEX)

install: bemcp-install
	$(COPY) fieldtrip/src/read_24bit.$(SUF) fieldtrip/fileio/private/
	$(COPY) fieldtrip/src/mxSerialize.$(SUF) fieldtrip/fileio/private/
	$(COPY) fieldtrip/src/mxDeserialize.$(SUF) fieldtrip/fileio/private/
	$(COPY) fieldtrip/src/read_16bit.$(SUF) fieldtrip/fileio/private/
	
	$(COPY) fieldtrip/src/read_ctf_shm.$(SUF) fieldtrip/fileio/private/ # Unix only
	$(COPY) fieldtrip/src/write_ctf_shm.$(SUF) fieldtrip/fileio/private/ # Unix only
	
	$(COPY) fieldtrip/src/meg_leadfield1.$(SUF) fieldtrip/forward/private/
	$(COPY) fieldtrip/src/ptriproj.$(SUF) fieldtrip/forward/private/
	$(COPY) fieldtrip/src/lmoutr.$(SUF) fieldtrip/forward/private/
	$(COPY) fieldtrip/src/plgndr.$(SUF) fieldtrip/forward/private/
	$(COPY) fieldtrip/src/solid_angle.$(SUF) fieldtrip/forward/private/
	$(COPY) fieldtrip/src/routlm.$(SUF) fieldtrip/forward/private/
	
	$(COPY) fieldtrip/src/solid_angle.$(SUF) fieldtrip/inverse/private/
	
	$(COPY) fieldtrip/src/ptriproj.$(SUF) fieldtrip/private/
	$(COPY) fieldtrip/src/lmoutr.$(SUF) fieldtrip/private/
	$(COPY) fieldtrip/src/plgndr.$(SUF) fieldtrip/private/
	$(COPY) fieldtrip/src/solid_angle.$(SUF) fieldtrip/private/

tarball: all bemcp-tarball
	$(TAR) cf spm_mex.tar $(FIELDTRIPMEX)
	$(TAR) Af spm_mex.tar $(BEMCP)/spm_mex.tar

bemcp bemcp-all bemcp-clean bemcp-distclean bemcp-install bemcp-tarball:
	@ $(MAKE) -C $(BEMCP) --no-print-directory $(subst -,,$(subst bemcp,,$@));

fieldtrip/src/%.$(SUF) : fieldtrip/src/%.c
	$(MEX) $< -outdir $(shell dirname $<) $(MEXEND)

fieldtrip/fileio/@uint64/%.$(SUF) : fieldtrip/fileio/@uint64/%.c
	$(MEX) $< -outdir $(shell dirname $<) $(MEXEND)

fieldtrip/src/lmoutr.$(SUF): fieldtrip/src/lmoutr.c fieldtrip/src/geometry.c fieldtrip/src/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/src/geometry.c $(MEXEND)

fieldtrip/src/routlm.$(SUF): fieldtrip/src/routlm.c fieldtrip/src/geometry.c fieldtrip/src/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/src/geometry.c $(MEXEND)

fieldtrip/src/ptriproj.$(SUF): fieldtrip/src/ptriproj.c fieldtrip/src/geometry.c fieldtrip/src/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/src/geometry.c $(MEXEND)

fieldtrip/src/solid_angle.$(SUF): fieldtrip/src/solid_angle.c fieldtrip/src/geometry.c fieldtrip/src/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/src/geometry.c $(MEXEND)

fieldtrip/src/ltrisect.$(SUF): fieldtrip/src/ltrisect.c fieldtrip/src/geometry.c fieldtrip/src/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/src/geometry.c $(MEXEND)

fieldtrip/src/plinproj.$(SUF): fieldtrip/src/plinproj.c fieldtrip/src/geometry.c fieldtrip/src/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/src/geometry.c $(MEXEND)
