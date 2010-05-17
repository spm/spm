#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 3890 2010-05-17 11:17:19Z guillaume $

include ../src/Makefile.var

FIELDTRIPMEX =\
	fieldtrip/fileio/private/mxDeserialize.$(SUF) fieldtrip/fileio/private/mxSerialize.$(SUF) \
	fieldtrip/fileio/private/read_16bit.$(SUF) fieldtrip/fileio/private/read_24bit.$(SUF) \
	fieldtrip/fileio/@uint64/max.$(SUF) fieldtrip/fileio/@uint64/min.$(SUF) \
	fieldtrip/fileio/@uint64/minus.$(SUF) fieldtrip/fileio/@uint64/plus.$(SUF) \
	fieldtrip/fileio/@uint64/rdivide.$(SUF) fieldtrip/fileio/@uint64/times.$(SUF) \
	fieldtrip/fileio/@uint64/abs.$(SUF) \
	\
	fieldtrip/forward/private/lmoutr.$(SUF) fieldtrip/forward/private/meg_leadfield1.$(SUF) \
	fieldtrip/forward/private/plgndr.$(SUF) fieldtrip/forward/private/ptriproj.$(SUF) \
	fieldtrip/forward/private/solid_angle.$(SUF) \
	\
	fieldtrip/inverse/private/solid_angle.$(SUF) \
	\
	fieldtrip/private/lmoutr.$(SUF) fieldtrip/private/ltrisect.$(SUF) \
	fieldtrip/private/plgndr.$(SUF) fieldtrip/private/plinproj.$(SUF) \
	fieldtrip/private/ptriproj.$(SUF) fieldtrip/private/routlm.$(SUF) \
	fieldtrip/private/solid_angle.$(SUF) fieldtrip/private/splint_gh.$(SUF) \
	\
	fieldtrip/fileio/private/read_ctf_shm.$(SUF) fieldtrip/fileio/private/write_ctf_shm.$(SUF) # Unix only

BEMCP = bemcp

all: $(FIELDTRIPMEX) $(FILEIOMEX) $(FORWINVMEX) bemcp-all

clean: bemcp-clean

distclean: clean bemcp-distclean
	$(DEL) $(FIELDTRIPMEX) $(FILEIOMEX) $(FORWINVMEX)

install: bemcp-install

tarball: all bemcp-tarball
	$(TAR) cf spm_mex.tar $(FIELDTRIPMEX) $(FILEIOMEX) $(FORWINVMEX)
	$(TAR) Af spm_mex.tar $(BEMCP)/spm_mex.tar

bemcp bemcp-all bemcp-clean bemcp-distclean bemcp-install bemcp-tarball:
	@ $(MAKE) -C $(BEMCP) --no-print-directory $(subst -,,$(subst bemcp,,$@));

%.$(SUF) : %.c
	$(MEX) $< -outdir $(shell dirname $<) $(MEXEND)

fieldtrip/forward/private/lmoutr.$(SUF): fieldtrip/forward/private/lmoutr.c fieldtrip/forward/private/geometry.c fieldtrip/forward/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/forward/private/geometry.c $(MEXEND)

fieldtrip/forward/private/ltrisect.$(SUF): fieldtrip/forward/private/ltrisect.c fieldtrip/forward/private/geometry.c fieldtrip/forward/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/forward/private/geometry.c $(MEXEND)

fieldtrip/forward/private/plinproj.$(SUF): fieldtrip/forward/private/plinproj.c fieldtrip/forward/private/geometry.c fieldtrip/forward/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/forward/private/geometry.c $(MEXEND)

fieldtrip/forward/private/ptriproj.$(SUF): fieldtrip/forward/private/ptriproj.c fieldtrip/forward/private/geometry.c fieldtrip/forward/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/forward/private/geometry.c $(MEXEND)

fieldtrip/forward/private/routlm.$(SUF): fieldtrip/forward/private/routlm.c fieldtrip/forward/private/geometry.c fieldtrip/forward/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/forward/private/geometry.c $(MEXEND)

fieldtrip/forward/private/solid_angle.$(SUF): fieldtrip/forward/private/solid_angle.c fieldtrip/forward/private/geometry.c fieldtrip/forward/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/forward/private/geometry.c $(MEXEND)


fieldtrip/private/lmoutr.$(SUF): fieldtrip/private/lmoutr.c fieldtrip/private/geometry.c fieldtrip/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/private/geometry.c $(MEXEND)

fieldtrip/private/ltrisect.$(SUF): fieldtrip/private/ltrisect.c fieldtrip/private/geometry.c fieldtrip/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/private/geometry.c $(MEXEND)

fieldtrip/private/plinproj.$(SUF): fieldtrip/private/plinproj.c fieldtrip/private/geometry.c fieldtrip/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/private/geometry.c $(MEXEND)

fieldtrip/private/ptriproj.$(SUF): fieldtrip/private/ptriproj.c fieldtrip/private/geometry.c fieldtrip/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/private/geometry.c $(MEXEND)

fieldtrip/private/routlm.$(SUF): fieldtrip/private/routlm.c fieldtrip/private/geometry.c fieldtrip/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/private/geometry.c $(MEXEND)

fieldtrip/private/solid_angle.$(SUF): fieldtrip/private/solid_angle.c fieldtrip/private/geometry.c fieldtrip/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) fieldtrip/private/geometry.c $(MEXEND)
