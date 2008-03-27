#!/usr/bin/env make -f
# General Makefile to compile SPM C-MEX files
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 1256 2008-03-27 20:42:31Z guillaume $
#
###############################################################################
#
# This Makefile has been tested under Linux, Windows, MacOS and SunOS.
# 
# If you have to tweak this Makefile to compile the SPM mex-files for your 
# platform, please send the details to <spm@fil.ion.ucl.ac.uk> so they can be 
# included here. 
#
# You can specify a particular platform with the following syntax:
#    make PLATFORM=Your_Platform
# The standard targets are 'all', 'clean', 'distclean' and 'install'.
#
# For a list of compatible compilers, see
#    http://www.mathworks.com/support/compilers/current_release/
#
# You can find some more help online on the SPM wikibook:
#  * Linux:
#      http://en.wikibooks.org/wiki/SPM/Installation_on_Linux
#      http://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Linux 
#  * Windows:
#      http://en.wikibooks.org/wiki/SPM/Installation_on_Windows
#      MinGW: http://www.mingw.org/
#      Gnumex: http://gnumex.sourceforge.net/
#  * SunOS:
#      http://en.wikibooks.org/wiki/SPM/Installation_on_SunOS
#  * MacOS:
#      http://en.wikibooks.org/wiki/SPM/Installation_on_Mac_OS
#      http://en.wikibooks.org/wiki/SPM/Installation_on_Mac_OS_(Intel)
#
###############################################################################

SHELL         = /bin/sh
MEXBIN        = mex
MEXOPTS       = -O
MEXEND        = 
CC            = cc
CFLAGS        = 
MOSUF         = o # mex output object suffix
CHMODIT       = chmod 644
AR            = ar rcs
COPY          = cp -f
DEL           = rm -f
MOVE          = mv -f

ifndef SUF
  ifndef PLATFORM
    PLATFORM  = $(shell uname)
  endif
  ifeq (SunOS,$(PLATFORM))
    SUF       = mexsol
    CC        = gcc
    CFLAGS    = -O3 -funroll-loops -DBIGENDIAN -fPIC
    MEXOPTS   = COPTIMFLAGS='-O3 -funroll-loops -DBIGENDIAN'
  endif
  ifeq (SunOS64,$(PLATFORM))
    SUF       = mexs64
    CC        = gcc
    CFLAGS    = -m64 -O3 -funroll-loops -DBIGENDIAN -fPIC
    MEXOPTS   = COPTIMFLAGS='-O3 -funroll-loops -DBIGENDIAN'
  endif
  ifeq (Linux,$(PLATFORM))
    HARDWARE  = $(shell uname -m)
    ifeq (i386,$(HARDWARE))
      SUF     = mexglx
      CC      = gcc
      CFLAGS  = -O3 -funroll-loops -fPIC -fexceptions
      MEXOPTS = COPTIMFLAGS='-O3 -funroll-loops -fexceptions'
    endif
    ifeq (i686,$(HARDWARE))
      SUF     = mexglx
      CC      = gcc
      CFLAGS  = -O3 -funroll-loops -fPIC -fexceptions
      MEXOPTS = COPTIMFLAGS='-O3 -funroll-loops -fexceptions'
    endif
    ifeq (x86_64,$(HARDWARE))
      SUF     = mexa64
      CC      = gcc
      CFLAGS  = -O3 -funroll-loops -fPIC -march=x86-64 -mfpmath=sse
      MEXOPTS = COPTIMFLAGS='-O3 -funroll-loops -fPIC -march=x86-64 -mfpmath=sse'
    endif
    ifndef SUF
      $(error Unknowm platform $(PLATFORM)-$(HARDWARE))
    endif
  endif
  ifeq (Darwin,$(PLATFORM))
    HARDWARE  = $(shell uname -p)
    ifeq (i386,$(HARDWARE))
       SUF       = mexmaci
       CFLAGS    = 
       MEXOPTS   = -O
    endif
    ifeq (powerpc,$(HARDWARE))
       SUF       = mexmac
       CFLAGS    = -DBIGENDIAN
       MEXOPTS   = -O -DBIGENDIAN
    endif
  endif
  ifeq (windows,$(PLATFORM))
    override PLATFORM = MINGW32
  endif
  ifeq (MINGW32,$(word 1,$(subst _, ,$(PLATFORM))))
    SUF       = mexw32
    CC        = gcc
    CFLAGS    = -O3 -funroll-loops -fomit-frame-pointer -DSPM_WIN32
    MEXBIN    = cmd /c "mex.bat
    MEXOPTS   = -O -DSPM_WIN32
    MEXEND    = "
    MOSUF     = obj
    CHMODIT   = chmod 755
  endif
  ifndef SUF
    $(error Unknowm platform $(PLATFORM))
  endif
endif

MEX           = $(MEXBIN) $(MEXOPTS)

###############################################################################
# Objects to go in the archive and mexfiles
###############################################################################

OBS    =\
	utils_uchar.$(SUF).o utils_short.$(SUF).o utils_int.$(SUF).o \
	utils_schar.$(SUF).o utils_ushort.$(SUF).o utils_uint.$(SUF).o\
	utils_float.$(SUF).o utils_double.$(SUF).o\
	utils_short_s.$(SUF).o utils_int_s.$(SUF).o\
	utils_ushort_s.$(SUF).o utils_uint_s.$(SUF).o\
	utils_float_s.$(SUF).o utils_double_s.$(SUF).o\
	spm_make_lookup.$(SUF).o spm_getdata.$(SUF).o spm_vol_access.$(SUF).o\
	spm_mapping.$(SUF).o

SPMMEX =\
	spm_sample_vol.$(SUF) spm_slice_vol.$(SUF) spm_brainwarp.$(SUF)\
	spm_add.$(SUF) spm_conv_vol.$(SUF) spm_render_vol.$(SUF)\
	spm_global.$(SUF) spm_resels_vol.$(SUF)\
	spm_bsplinc.$(SUF) spm_bsplins.$(SUF) spm_bias_mex.$(SUF)\
	spm_atranspa.$(SUF) spm_unlink.$(SUF) spm_hist.$(SUF)\
	spm_krutil.$(SUF) spm_project.$(SUF) spm_hist2.$(SUF)\
	spm_dilate_erode.$(SUF) spm_bwlabel.$(SUF) spm_get_lm.$(SUF)\
	spm_digamma.$(SUF) mat2file.$(SUF) file2mat.$(SUF)\
	spm_invdef.$(SUF)

###############################################################################
# The main ways to run make
###############################################################################

all: verb.$(SUF) $(SPMMEX)
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        FINISHED"
	@ echo "_____________________________________________________________"
	@ echo ""

clean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting object (.o) files"
	@ echo "_____________________________________________________________"
	@ echo ""
	$(DEL) $(OBS)

distclean: clean
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting mex and archive (.a) files"
	@ echo "_____________________________________________________________"
	@ echo ""
	$(DEL) $(SPMMEX) spm_vol_utils.$(SUF).a

archive: spm_vol_utils.$(SUF).a

install: $(SPMMEX)
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Installing mex files"
	@ echo "_____________________________________________________________"
	@ echo ""	
	$(COPY) $(SPMMEX) ..

###############################################################################
# Compile spm_vol_utils.c with various flags
###############################################################################

spm_vol_utils.$(SUF).a: $(OBS)
	$(DEL) $@
	$(AR) $@ $(OBS)
	@ $(CHMODIT) $@

UTILS = spm_vol_utils.c spm_make_lookup.h spm_getdata.h

utils_uchar.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_UNSIGNED_CHAR
	@ $(CHMODIT) $@

utils_short.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_SIGNED_SHORT
	@ $(CHMODIT) $@

utils_int.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_SIGNED_INT
	@ $(CHMODIT) $@

utils_schar.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_SIGNED_CHAR
	@ $(CHMODIT) $@

utils_ushort.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_UNSIGNED_SHORT
	@ $(CHMODIT) $@

utils_uint.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_UNSIGNED_INT
	@ $(CHMODIT) $@

utils_float.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_FLOAT
	@ $(CHMODIT) $@

utils_double.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_DOUBLE
	@ $(CHMODIT) $@

utils_short_s.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_int_s.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_ushort_s.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_uint_s.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_float_s.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_double_s.$(SUF).o: $(UTILS)
	$(CC) -c $(CFLAGS) -o $@ spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP
	@ $(CHMODIT) $@

###############################################################################
# Compile a few additional C routines for linking
###############################################################################

%.$(SUF).o : %.c
	$(CC) -c $(CFLAGS) -o $@ $<
	@ $(CHMODIT) $@

spm_getdata.$(SUF).o: spm_getdata.c spm_getdata.h
	$(CC) -c $(CFLAGS) -o $@ spm_getdata.c
	@ $(CHMODIT) $@

spm_vol_access.$(SUF).o: spm_vol_access.c spm_vol_access.h spm_datatypes.h
	$(CC) -c $(CFLAGS) -o $@ spm_vol_access.c
	@ $(CHMODIT) $@

spm_make_lookup.$(SUF).o: spm_make_lookup.c spm_make_lookup.h
	$(CC) -c $(CFLAGS) -o $@ spm_make_lookup.c
	@ $(CHMODIT) $@

spm_mapping.$(SUF).o: spm_mapping.c spm_mapping.h spm_vol_access.h spm_datatypes.h
	$(MEX) -c spm_mapping.c $(MEXEND)
	$(MOVE) spm_mapping.$(MOSUF) $@
	@ $(CHMODIT) $@

###############################################################################
# Compile the mex files themselves
###############################################################################

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
	@ $(CHMODIT) $@

spm_add.$(SUF): spm_add.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_add.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_bias_mex.$(SUF): spm_bias_mex.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_datatypes.h
	$(MEX) spm_bias_mex.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_brainwarp.$(SUF): spm_brainwarp.c spm_matfuns.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_brainwarp.c spm_vol_utils.$(SUF).a spm_matfuns.c $(MEXEND)
	@ $(CHMODIT) $@

spm_bsplinc.$(SUF): spm_bsplinc.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_bsplinc.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_conv_vol.$(SUF): spm_conv_vol.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h spm_datatypes.h
	$(MEX) spm_conv_vol.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_global.$(SUF): spm_global.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_global.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_render_vol.$(SUF): spm_render_vol.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_render_vol.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_resels_vol.$(SUF): spm_resels_vol.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_resels_vol.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_sample_vol.$(SUF): spm_sample_vol.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_sample_vol.c spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

spm_slice_vol.$(SUF): spm_slice_vol.c  spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_slice_vol.c  spm_vol_utils.$(SUF).a $(MEXEND)
	@ $(CHMODIT) $@

file2mat.$(SUF): file2mat.c
	$(MEX) file2mat.c $(MEXEND)
	@ $(CHMODIT) $@

mat2file.$(SUF): mat2file.c
	$(MEX) mat2file.c $(MEXEND)
	@ $(CHMODIT) $@

spm_atranspa.$(SUF): spm_atranspa.c 
	$(MEX) spm_atranspa.c $(MEXEND)
	@ $(CHMODIT) $@

spm_bsplins.$(SUF): spm_bsplins.c
	$(MEX) spm_bsplins.c $(MEXEND)
	@ $(CHMODIT) $@
	
spm_bwlabel.$(SUF): spm_bwlabel.c
	$(MEX) spm_bwlabel.c $(MEXEND)
	@ $(CHMODIT) $@

spm_digamma.$(SUF): spm_digamma.c
	$(MEX) spm_digamma.c $(MEXEND)
	@ $(CHMODIT) $@

spm_dilate_erode.$(SUF): spm_dilate_erode.c
	$(MEX) spm_dilate_erode.c $(MEXEND)
	@ $(CHMODIT) $@

spm_get_lm.$(SUF): spm_get_lm.c
	$(MEX) spm_get_lm.c $(MEXEND)
	@ $(CHMODIT) $@

spm_hist.$(SUF): spm_hist.c
	$(MEX) spm_hist.c $(MEXEND)
	@ $(CHMODIT) $@

spm_hist2.$(SUF): spm_hist2.c
	$(MEX) spm_hist2.c $(MEXEND)
	@ $(CHMODIT) $@

spm_invdef.$(SUF): spm_invdef.c
	$(MEX) spm_invdef.c $(MEXEND)
	@ $(CHMODIT) $@

spm_krutil.$(SUF): spm_krutil.c
	$(MEX) spm_krutil.c $(MEXEND)
	@ $(CHMODIT) $@

spm_project.$(SUF): spm_project.c
	$(MEX) spm_project.c $(MEXEND)
	@ $(CHMODIT) $@	

spm_unlink.$(SUF): spm_unlink.c
	$(MEX) spm_unlink.c $(MEXEND)
	@ $(CHMODIT) $@

###############################################################################
# Assorted architecture dependent messages
###############################################################################

verb.mexw32:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "     Windows compilation (32 bit) - using $(CC)/`uname`"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexw64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "     Windows compilation (64 bit) - using $(CC)/`uname`"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexglx:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Linux compilation (x86-32) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexa64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Linux compilation (x86-64) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmac:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (MacOS X, PowerPC) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmaci:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (MacOS X, Intel) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexsol:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (Solaris 32 bit) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexs64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (Solaris 64 bit) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""
