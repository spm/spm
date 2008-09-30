#!/usr/bin/env make -f
# General Makefile to compile SPM C-MEX files
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 2267 2008-09-30 21:13:55Z guillaume $
#
###############################################################################
#
# This Makefile has been tested under Linux, Windows, MacOS and SunOS.
# 
# If you have to tweak this Makefile or Makefile.var to compile the SPM 
# mex-files for your platform, please send the details to <spm@fil.ion.ucl.ac.uk> 
# so they can be included here. 
#
# To compile and install SPM, type the following in a Unix console:
# >  make && make install
#
# You can specify a particular platform with the following syntax:
# >  make PLATFORM=Your_Platform
# The standard targets are 'all', 'clean', 'distclean' and 'install'.
#
# For a list of compatible compilers, see
#    http://www.mathworks.com/support/compilers/current_release/
#
###############################################################################

include Makefile.var

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
	spm_unlink.$(SUF) spm_existfile.$(SUF) spm_hist.$(SUF)\
	spm_krutil.$(SUF) spm_project.$(SUF) spm_hist2.$(SUF)\
	spm_dilate_erode.$(SUF) spm_bwlabel.$(SUF) spm_get_lm.$(SUF)\
	spm_invdef.$(SUF) mat2file.$(SUF) file2mat.$(SUF) spm_voronoi.$(SUF)

SUBDIRS=\
	toolbox/DARTEL \
	toolbox/FieldMap \
	toolbox/HDW \
	@xmltree/private

###############################################################################
# The main ways to run make
###############################################################################

all: verb.$(SUF) $(SPMMEX) verb.end

clean: verb.clean
	$(DEL) $(OBS)

distclean: clean verb.distclean
	$(DEL) $(SPMMEX) spm_vol_utils.$(SUF).a

archive: spm_vol_utils.$(SUF).a

install: verb.install $(SPMMEX)
	$(COPY) $(SPMMEX) ..
	$(MOVE) ../file2mat.$(SUF) ../mat2file.$(SUF) ../@file_array/private

toolbox toolbox-all toolbox-clean toolbox-distclean toolbox-install toolbox-tarball:
	@@for p in $(SUBDIRS); do \
	echo "_____________________________________________________________";\
	echo "";\
	echo "          In " ../$$p; \
	echo "_____________________________________________________________";\
	echo ""	;\
	$(MAKE) -C ../$$p --no-print-directory $(subst -,,$(subst toolbox,,$@));\
	done

tarball: all toolbox-all toolbox-tarball verb.archive
	$(TAR) cf spm_mex.tar $(SPMMEX)
	@@for p in $(SUBDIRS); do \
		$(TAR) Af spm_mex.tar  ../$$p/spm_mex.tar;\
	done
	$(ZIP) spm_mex.tar

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

spm_bsplins.$(SUF): spm_bsplins.c
	$(MEX) spm_bsplins.c $(MEXEND)
	@ $(CHMODIT) $@
	
spm_bwlabel.$(SUF): spm_bwlabel.c
	$(MEX) spm_bwlabel.c $(MEXEND)
	@ $(CHMODIT) $@

spm_dilate_erode.$(SUF): spm_dilate_erode.c
	$(MEX) spm_dilate_erode.c $(MEXEND)
	@ $(CHMODIT) $@

spm_existfile.$(SUF): spm_existfile.c
	$(MEX) spm_existfile.c $(MEXEND)
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

spm_voronoi.$(SUF): spm_voronoi.c
	$(MEX) spm_voronoi.c $(MEXEND)
	@ $(CHMODIT) $@

###############################################################################
# General messages
###############################################################################

verb.end:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        FINISHED"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.clean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting object (.o) files"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.distclean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting mex and archive (.a) files"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.install:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Installing mex files"
	@ echo "_____________________________________________________________"
	@ echo ""	

verb.archive:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Creating archive spm_mex.tar.gz"
	@ echo "_____________________________________________________________"
	@ echo ""	

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
	@ echo "        Unix compilation (MacOS X, Intel 32 bit) - using $(CC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmaci64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (MacOS X, Intel 64 bit) - using $(CC)"
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
