#!/usr/bin/env make -f
# Manual Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2016-2018 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 7485 2018-11-14 10:31:57Z guillaume $

include ../src/Makefile.var

PDFLATEX = pdflatex -interaction=batchmode
BIBTEX   = bibtex
ifeq (octave,$(PLATFORM))
  SPMEXE = ../bin/spm12-octave
else
  SPMEXE = ../bin/spm12-matlab
endif

SPMPDF   = manual.pdf ReleaseNotes.pdf

all: $(SPMPDF)
	@:

clean:
	$(DEL) *.toc *.aux */*.aux *.log *.out *.bbl *.blg *.maf *.mtc*
        
distclean: clean
	$(DEL) $(SPMPDF)

install:
	@:

%.pdf: %.tex
	$(call makepdf, $*)

manual.pdf: manual.tex spm_manual.tex
	$(call makepdf, "manual")

spm_manual.tex:
	$(call verb, "Creating TEX files")
	@ ${SPMEXE} eval "spm_jobman initcfg;spm_make_manual;"
	$(call verb, "TEX files: done")

define makepdf
	$(call verb, "Building $(1).pdf")
	@ $(PDFLATEX) $(1) && $(BIBTEX) $(1) && $(PDFLATEX) $(1) && $(PDFLATEX) $(1)
        $(call verb, "$(1).pdf: done")
endef
