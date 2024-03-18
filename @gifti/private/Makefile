#!/usr/bin/env make -f
# GIfTI Makefile
#
# Copyright (C) 2008-2023 Wellcome Centre for Human Neuroimaging

include ../../src/Makefile.var

SPMMEX = zstream.$(MEXEXT) xml_parser.$(MEXEXT) base64.$(MEXEXT)

ifeq (mex,$(MEXEXT))
  export CFLAGS = $(shell $(MEX) -p CFLAGS) -std=c99
else
  ifeq (windows,$(PLATFORM))
    MEXOPTS += CFLAGS=\"$$CFLAGS -std=c99\"
  else
    MEXOPTS += CFLAGS='$$CFLAGS -std=c99'
  endif
endif

all: $(SPMMEX)
	@:

clean:
	@:
        
distclean: clean
	$(DEL) $(SPMMEX)

install:
	@:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(MEXEXT) : %.c
	$(MEX) $< $(MEXEND)

xml_parser.$(MEXEXT): xml_parser.c yxml.c yxml.h
	$(MEX) $< yxml.c $(MEXEND)
