#!/usr/bin/env make -f
# Makefile called by {SPM}/src/Makefile to compile SPM external's
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 2990 2009-03-28 17:33:52Z guillaume $

include ../src/Makefile.var

FIELDTRIPMEX =\
	fieldtrip/@config/private/deepcopy.$(SUF) fieldtrip/@config/private/increment.$(SUF) \
	fieldtrip/@config/private/reset.$(SUF) \
	fieldtrip/@uint64/max.$(SUF) fieldtrip/@uint64/min.$(SUF) \
	fieldtrip/@uint64/minus.$(SUF) fieldtrip/@uint64/plus.$(SUF) \
	fieldtrip/@uint64/rdivide.$(SUF) fieldtrip/@uint64/times.$(SUF) \
	fieldtrip/@uint64/abs.$(SUF) \
	fieldtrip/private/lmoutr.$(SUF) fieldtrip/private/ltrisect.$(SUF) \
	fieldtrip/private/meg_leadfield1.$(SUF) fieldtrip/private/plgndr.$(SUF) \
	fieldtrip/private/plinproj.$(SUF) fieldtrip/private/ptriproj.$(SUF) \
	fieldtrip/private/read_24bit.$(SUF) fieldtrip/private/routlm.$(SUF) \
	fieldtrip/private/solid_angle.$(SUF) fieldtrip/private/splint_gh.$(SUF) 

FILEIOMEX =\
	fileio/@uint64/max.$(SUF) fileio/@uint64/min.$(SUF) \
	fileio/@uint64/minus.$(SUF) fileio/@uint64/plus.$(SUF) \
	fileio/@uint64/rdivide.$(SUF) fileio/@uint64/times.$(SUF) \
	fileio/@uint64/abs.$(SUF) \
	fileio/private/read_24bit.$(SUF) \
	fileio/private/mxSerialize.$(SUF) fileio/private/mxDeserialize.$(SUF) \
	fileio/private/read_ctf_shm.$(SUF) fileio/private/write_ctf_shm.$(SUF) # Unix only

FORWINVMEX =\
	forwinv/private/lmoutr.$(SUF) forwinv/private/ltrisect.$(SUF) \
	forwinv/private/meg_leadfield1.$(SUF) forwinv/private/plgndr.$(SUF) \
	forwinv/private/plinproj.$(SUF) forwinv/private/ptriproj.$(SUF) \
	forwinv/private/routlm.$(SUF) forwinv/private/solid_angle.$(SUF)

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


forwinv/private/lmoutr.$(SUF): forwinv/private/lmoutr.c forwinv/private/geometry.c forwinv/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) forwinv/private/geometry.c $(MEXEND)

forwinv/private/ltrisect.$(SUF): forwinv/private/ltrisect.c forwinv/private/geometry.c forwinv/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) forwinv/private/geometry.c $(MEXEND)

forwinv/private/plinproj.$(SUF): forwinv/private/plinproj.c forwinv/private/geometry.c forwinv/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) forwinv/private/geometry.c $(MEXEND)

forwinv/private/ptriproj.$(SUF): forwinv/private/ptriproj.c forwinv/private/geometry.c forwinv/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) forwinv/private/geometry.c $(MEXEND)

forwinv/private/routlm.$(SUF): forwinv/private/routlm.c forwinv/private/geometry.c forwinv/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) forwinv/private/geometry.c $(MEXEND)

forwinv/private/solid_angle.$(SUF): forwinv/private/solid_angle.c forwinv/private/geometry.c forwinv/private/geometry.h
	$(MEX) $< -outdir $(shell dirname $<) forwinv/private/geometry.c $(MEXEND)
