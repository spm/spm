#!make -f
#
# %W% John Ashburner & Matthew Brett %E%
# $Id: Makefile,v 2.6 2002-08-08 17:54:13 john Exp $
#
###############################################################################
#
# Suggestions for how to make this file a bit more elegant are welcome.  So far
# it works under SunOS, Linux (at the FIL) and some Windows systems.
# $Log: not supported by cvs2svn $
#
###############################################################################

###############################################################################
# Default make
###############################################################################
unknown:
	@ make `uname`
tidy:
	@ make clean.`uname`

###############################################################################
# Architecture specific makes
###############################################################################

MEX = mex -O
CC  = cc
# mex output object suffix
MOSUF = o
CHMODIT = chmod 644 
ADDED_OBS=

SunOS:
	make all SUF=mexsol  CC="cc -xO5"                   MEX="mex COPTIMFLAGS=-xO5"
Linux:
	make all SUF=mexglx  CC="gcc -O3 -funroll-loops"    MEX="mex COPTIMFLAGS='-O3 -funroll-loops'"
HP-UX:
	make all SUF=mexhp7  CC="cc  -O +z -Ae +DAportable" MEX="mex COPTIMFLAGS=-O"
IRIX:
	make all SUF=mexsg   CC="cc  -O -mips2"             MEX="mex -O"
IRIX64:
	make all SUF=mexsg64 CC="cc  -O -mips4 -64"         MEX="mex -O"
AIX:
	make all SUF=mexrs6  MEX="mex -O"
OSF1:
	make all SUF=mexaxp  MEX="mex -O"
windows:
	make all SUF=dll     CC="gcc -mno-cygwin -DSPM_WIN32" MEX="mex.bat -DSPM_WIN32" MOSUF=obj CHMODIT="echo > null" ADDED_OBS=win32mmap.dll.o

###############################################################################
# Architecture specific cleaning
###############################################################################
clean.SunOS:
	make clean SUF=mexsol
clean.Linux:
	make clean SUF=mexglx
clean.HP-UX:
	make clean SUF=mexhp7
clean.IRIX:
	make clean SUF=mexsg
clean.IRIX64:
	make clean SUF=mexsg64
clean.AIX:
	make clean SUF=mexrs6
clean.OSF1:
	make clean SUF=mexaxp
clean.windows:
	make clean SUF=dll

###############################################################################
# Objects to go in the archive and mexfiles
###############################################################################

SUF = unknown

OBS =	utils_uchar.$(SUF).o utils_short.$(SUF).o utils_int.$(SUF).o \
	utils_schar.$(SUF).o utils_ushort.$(SUF).o utils_uint.$(SUF).o\
	utils_float.$(SUF).o utils_double.$(SUF).o\
	utils_short_s.$(SUF).o utils_int_s.$(SUF).o\
	utils_ushort_s.$(SUF).o utils_uint_s.$(SUF).o\
	utils_float_s.$(SUF).o utils_double_s.$(SUF).o\
	spm_make_lookup.$(SUF).o spm_getdata.$(SUF).o spm_vol_access.$(SUF).o\
	spm_mapping.$(SUF).o $(ADDED_OBS)

SPMMEX =\
	spm_sample_vol.$(SUF) spm_slice_vol.$(SUF) spm_brainwarp.$(SUF)\
	spm_add.$(SUF) spm_conv_vol.$(SUF) spm_render_vol.$(SUF)\
	spm_global.$(SUF) spm_resels_vol.$(SUF) spm_getxyz.$(SUF)\
	spm_atranspa.$(SUF) spm_list_files.$(SUF) spm_unlink.$(SUF)\
	spm_krutil.$(SUF) spm_project.$(SUF) spm_hist2.$(SUF) spm_max.$(SUF)\
	spm_clusters.$(SUF) spm_bsplinc.$(SUF) spm_bsplins.$(SUF)\
	spm_bias_mex.$(SUF)

###############################################################################
# The main ways to run make
###############################################################################

all: verb.$(SUF) $(SPMMEX)
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        FINISHED"
	@ echo "_____________________________________________________________"
	@ echo ""

very_clean: clean
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting mex and archive (.a) files"
	@ echo "_____________________________________________________________"
	@ echo ""
	rm -f $(SPMMEX) spm_vol_utils.$(SUF).a

clean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting object (.o) files"
	@ echo "_____________________________________________________________"
	@ echo ""
	rm -f $(OBS)

archive: spm_vol_utils.$(SUF).a

###############################################################################
# Compile spm_vol_utils.c with various flags
###############################################################################
CC  = cc                                                                                                              
spm_vol_utils.$(SUF).a: $(OBS)
	rm -f $@
	ar rcv $@ $(OBS)
	@ $(CHMODIT) $@

UTILS=spm_vol_utils.c spm_sys_deps.h spm_make_lookup.h spm_getdata.h

utils_uchar.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_UNSIGNED_CHAR
	@ $(CHMODIT) $@

utils_short.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_SIGNED_SHORT
	@ $(CHMODIT) $@

utils_int.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_SIGNED_INT
	@ $(CHMODIT) $@

utils_schar.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_SIGNED_CHAR
	@ $(CHMODIT) $@

utils_ushort.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_UNSIGNED_SHORT
	@ $(CHMODIT) $@

utils_uint.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_UNSIGNED_INT
	@ $(CHMODIT) $@

utils_float.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_FLOAT
	@ $(CHMODIT) $@

utils_double.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_DOUBLE
	@ $(CHMODIT) $@

utils_short_s.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_int_s.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_ushort_s.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_uint_s.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_float_s.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP
	@ $(CHMODIT) $@

utils_double_s.$(SUF).o: $(UTILS)
	$(CC) -c -o $@ spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP
	@ $(CHMODIT) $@

###############################################################################
# Compile a few additional C routines for linking
###############################################################################

%.$(SUF).o : %.c spm_sys_deps.h
	$(CC) -c -o $@ $<
	@ $(CHMODIT) $@

spm_vol_access.$(SUF).o:  spm_vol_access.c spm_vol_access.h spm_datatypes.h

spm_make_lookup.$(SUF).o: spm_make_lookup.c spm_sys_deps.h

spm_mapping.$(SUF).o:     spm_mapping.c spm_sys_deps.h spm_mapping.h spm_datatypes.h
	$(MEX) -c spm_mapping.c
	mv spm_mapping.$(MOSUF) $@
	@ $(CHMODIT) $@

win32mmap.$(SUF).o: win32mmap.c win32mmap.h
	$(CC) -c -o $@ win32mmap.c
	@ $(CHMODIT) $@

###############################################################################
# Compile the mex files themselves
###############################################################################
MEX = mex -O
%.$(SUF) : %.c spm_sys_deps.h
	$(MEX) $<
	@ $(CHMODIT) $@

spm_add.$(SUF): spm_add.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_add.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_brainwarp.$(SUF): spm_brainwarp.c  spm_vol_utils.$(SUF).a spm_matfuns.c\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_brainwarp.c spm_vol_utils.$(SUF).a spm_matfuns.c
	@ $(CHMODIT) $@

spm_bsplinc.$(SUF): spm_bsplinc.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h spm_datatypes.h
	$(MEX) spm_bsplinc.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_bsplins.$(SUF): spm_bsplins.c spm_sys_deps.h
	$(MEX) spm_bsplins.c

spm_conv_vol.$(SUF): spm_conv_vol.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h spm_datatypes.h
	$(MEX) spm_conv_vol.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_getxyz.$(SUF): spm_getxyz.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_getxyz.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_global.$(SUF): spm_global.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_global.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_hist2.$(SUF):      spm_hist2.c spm_sys_deps.h

spm_krutil.$(SUF): spm_krutil.c spm_sys_deps.h

spm_list_files.$(SUF): spm_list_files.c spm_sys_deps.h

spm_project.$(SUF): spm_project.c spm_sys_deps.h

spm_render_vol.$(SUF): spm_render_vol.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_render_vol.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_resels_vol.$(SUF): spm_resels_vol.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_resels_vol.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_sample_vol.$(SUF): spm_sample_vol.c spm_vol_utils.$(SUF).a spm_mapping.h
	$(MEX) spm_sample_vol.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_slice_vol.$(SUF): spm_slice_vol.c  spm_vol_utils.$(SUF).a spm_mapping.h
	$(MEX) spm_slice_vol.c  spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_bias_mex.$(SUF): spm_bias_mex.c spm_vol_utils.$(SUF).a spm_mapping.h
	$(MEX) spm_bias_mex.c spm_vol_utils.$(SUF).a -DIGNORE_ZEROS
	@ $(CHMODIT) $@

###############################################################################
# Assorted architecture dependent messages
###############################################################################

verb.unknown:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "PROBLEM: Dont know how to do this."
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexhp7:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "unix compile for hpux cc, and maybe aix cc"
	@ echo ""
	@ echo "Under HPUX 10.20 with MATLAB 5.2.1 and gcc, you may wish"
	@ echo "to modify this Makefile to say something like:"
	@ echo '  CC=gcc -O -fpic'
	@ echo '  MEX=mex COPTIMFLAGS=-O -f gccopts.sh'
	@ echo "where the gccopts.sh file is modified to remove the +z"
	@ echo "(which is used with version 9.0 and possibly also 10.0)."
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexsg:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "Feedback from users with R10000 O2 and R10000 Indigo2 systems"
	@ echo "running IRIX6.5 suggests that the cmex program with Matlab 5.x"
	@ echo "compiles with the old 32bit (o32) instruction set (MIPS2) only,"
	@ echo "while cc by default compiles with the new32 bit (n32 or MIPS4)."
	@ echo "Matlab 5.x only likes o32 for O2 R10000 systems."
	@ echo ""
	@ echo "We also suggest you modify your options file mexopts.sh in"
	@ echo 'the sgi section: change LD="ld" to LD="ld -o32"'
	@ echo "this tells the linker to use o32 instead of n32."
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexsg64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "not optimised sgi 64 bit compile for CC"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.dll:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "Windows compile with gcc/mingw"
	@ echo "see http://www.mrc-cbu.cam.ac.uk/Imaging/gnumex20.html"
	@ echo "for instructions about installing gcc/mingw for"
	@ echo "compiling Mex files."
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexsol:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "Unix compile for Sun cc"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexglx:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "Linux compilation - using gcc"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexaxp:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "keep your fingers crossed"
	@ echo "_____________________________________________________________"
	@ echo ""

