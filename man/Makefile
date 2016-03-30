#!/usr/bin/env make -f
# Manual Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 6761 2016-03-30 15:11:20Z guillaume $

include ../src/Makefile.var

PDFLATEX = pdflatex -interaction=batchmode
BIBTEX   = bibtex

SPMPDF   = manual.pdf ReleaseNotes.pdf

all: verb.pdf $(SPMPDF) verb.pdf.end
	@:

clean:
	$(DEL) *.toc *.aux */*.aux *.log *.out *.bbl *.blg *.maf *.mtc*
        
distclean: clean
	$(DEL) $(SPMPDF)

install:
	@:

%.pdf: %.tex
	@ $(PDFLATEX) $* && $(BIBTEX) $* && $(PDFLATEX) $* && $(PDFLATEX) $*

verb.pdf:
	$(call verb, "Creating PDF files")

verb.pdf.end:
	$(call verb, "PDF files: done")
