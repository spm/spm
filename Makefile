#!make -f
#
# %W% John Ashburner %E%
# $Id$
#
###############################################################################
#
# Suggestions for how to make this file a bit more elegant are welcome.  So far
# it works under SunOS and Linux (at the FIL) and some Windows systems
#
# Thanks to Matthew, Darren, Alle Meije and others for various suggestions.
#
# For a list of compatible compilers, see
# http://www.mathworks.com/support/tech-notes/1600/1601.html
# http://www.mathworks.com/support/tech-notes/1600/1601_files/1601_65.html
#
# $Log$
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
RANLIB  = 
ADDED_OBS=
ADDED_MEX=

SunOS:
	make all SUF=mexsol  CC="cc -xO5 -DBIGENDIAN"\
	         MEX="mex COPTIMFLAGS=-xO5 -DBIGENDIAN"
SunOS.gcc:
# Assumes that gccopts.sh has been copied to the current directory
	make all SUF=mexsol  CC="gcc -O3 -funroll-loops -DBIGENDIAN -fPIC"\
	       MEX="mex COPTIMFLAGS='-O3 -funroll-loops -DBIGENDIAN' -f gccopts.sh"
Linux:
	make all SUF=mexglx  CC="gcc -O3 -funroll-loops -fPIC -fexceptions"\
	       MEX="mex COPTIMFLAGS='-O3 -funroll-loops -fexceptions'"
Linux.A64:
# The '-fPIC' option is nexessary to allow the linking proces to complete. 
# '-march=x86-64' provides generic optimisations for both Opteron and 64bit Xeon.
# If the code is running on a 64bit Xeon you can change '-march=x86-64' to 
# '-march=nocona', when running on an Opteron, change it to '-march=opteron'.
# Also possibly use '-march=k8' for 64 bit Athlon.
	make all SUF=mexa64  CC="gcc -O3 -funroll-loops -fPIC -march=x86-64 -mfpmath=sse"\
	       MEX="mex COPTIMFLAGS='-O3 -funroll-loops -fPIC -march=x86-64 -mfpmath=sse'"
HP-UX:
	make all SUF=mexhp7  CC="cc  -O +z -Ae +DAportable -DBIGENDIAN"
	       MEX="mex COPTIMFLAGS=-O -DBIGENDIAN"
IRIX:
	make all SUF=mexsg   CC="cc  -O -n32 -dont_warn_unused -OPT:IEEE_NaN_inf=ON -DBIGENDIAN"\
	       MEX="mex COPTIMFLAGS='-O' -DBIGENDIAN"
IRIX64:
	make all SUF=mexsg64   CC="cc  -O -mips4 -64 -DBIGENDIAN"\
	       MEX="mex COPTIMFLAGS='-O' -DBIGENDIAN"
AIX:
	make all SUF=mexrs6 CC="cc -DBIGENDIAN"\
	       MEX = "mex -O -DBIGENDIAN"
OSF1:
	make all SUF=mexaxp CC="cc -DBIGENDIAN"\
	       MEX = "mex -O -DBIGENDIAN"
MAC:
	make all SUF=mexmac RANLIB="ranlib spm_vol_utils.mexmac.a" CC="cc -DBIGENDIAN"\
	       MEX="mex -O -DBIGENDIAN"
#windows:
#	make all  SUF=dll MOSUF=obj CHMODIT="echo > null"  ADDED_OBS=win32mmap.dll.o ADDED_MEX=spm_win32utils.dll
windows:
# Consider adding either of the following, depending on your platform:
#       -march=pentium3
#       -march=pentium4
	make all SUF=dll CC="gcc -mno-cygwin -O3 -funroll-loops -fomit-frame-pointer -march=pentium4 -mfpmath=sse -DSPM_WIN32"\
	         MEX="mex.bat -DSPM_WIN32" MOSUF=obj CHMODIT="echo > NUL" ADDED_OBS=win32mmap.dll.o ADDED_MEX=spm_win32utils.dll

###############################################################################
# Architecture specific cleaning
###############################################################################
clean.SunOS:
	make clean SUF=mexsol
clean.SunOS.gcc:
	make clean SUF=mexsol
clean.Linux:
	make clean SUF=mexglx
clean.Linux.5:
	make clean SUF=mexlx
clean.Linux.A64:
	make clean SUF=mexa64
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
clean.MAC:
	make clean SUF=mexmac
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
	spm_global.$(SUF) spm_resels_vol.$(SUF)\
	spm_atranspa.$(SUF) spm_list_files.$(SUF) spm_unlink.$(SUF)\
	spm_krutil.$(SUF) spm_project.$(SUF) spm_hist2.$(SUF)\
	spm_bsplinc.$(SUF) spm_bsplins.$(SUF)\
	spm_bias_mex.$(SUF) spm_dilate.$(SUF)\
	spm_bwlabel.$(SUF) spm_get_lm.$(SUF)\
	spm_digamma.$(SUF) mat2file.$(SUF) file2mat.$(SUF)\
	$(ADDED_MEX)

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
                                                                                                          
spm_vol_utils.$(SUF).a: $(OBS)
	rm -f $@
	ar rcv $@ $(OBS)
	$(RANLIB)
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

spm_getdata.$(SUF).o: spm_getdata.c spm_sys_deps.h
	$(CC) -c -o $@ spm_getdata.c
	@ $(CHMODIT) $@

spm_vol_access.$(SUF).o:  spm_vol_access.c spm_vol_access.h spm_datatypes.h
	$(CC) -c -o $@ spm_vol_access.c
	@ $(CHMODIT) $@

spm_make_lookup.$(SUF).o: spm_make_lookup.c spm_sys_deps.h
	$(CC) -c -o $@ spm_make_lookup.c
	@ $(CHMODIT) $@

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
	@ $(CHMODIT) $@

spm_conv_vol.$(SUF): spm_conv_vol.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h spm_datatypes.h
	$(MEX) spm_conv_vol.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_global.$(SUF): spm_global.c spm_vol_utils.$(SUF).a\
		spm_sys_deps.h spm_mapping.h
	$(MEX) spm_global.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_atranspa.$(SUF): spm_atranspa.c spm_sys_deps.h
	$(MEX) spm_atranspa.c
	@ $(CHMODIT) $@

file2mat.$(SUF): file2mat.c
	$(MEX) file2mat.c
	@ $(CHMODIT) $@

mat2file.$(SUF): mat2file.c
	$(MEX) mat2file.c
	@ $(CHMODIT) $@

spm_digamma.$(SUF): spm_digamma.c
	$(MEX) spm_digamma.c
	@ $(CHMODIT) $@

spm_unlink.$(SUF): spm_unlink.c spm_sys_deps.h
	$(MEX) spm_unlink.c
	@ $(CHMODIT) $@

spm_hist2.$(SUF): spm_hist2.c spm_sys_deps.h
	$(MEX) spm_hist2.c
	@ $(CHMODIT) $@

spm_krutil.$(SUF): spm_krutil.c spm_sys_deps.h
	$(MEX) spm_krutil.c
	@ $(CHMODIT) $@

spm_dilate.$(SUF): spm_dilate.c
	$(MEX) spm_dilate.c
	@ $(CHMODIT) $@

spm_bwlabel.$(SUF): spm_bwlabel.c
	$(MEX) spm_bwlabel.c
	@ $(CHMODIT) $@

spm_get_lm.$(SUF): spm_get_lm.c
	$(MEX) spm_get_lm.c
	@ $(CHMODIT) $@

spm_list_files.$(SUF): spm_list_files.c spm_sys_deps.h
	$(MEX) spm_list_files.c
	@ $(CHMODIT) $@

spm_project.$(SUF): spm_project.c spm_sys_deps.h
	$(MEX) spm_project.c
	@ $(CHMODIT) $@

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
	$(MEX) spm_bias_mex.c spm_vol_utils.$(SUF).a
	@ $(CHMODIT) $@

spm_win32utils.$(SUF): spm_win32utils.c 
	$(MEX) spm_win32utils.c
	@ $(CHMODIT) $@

###############################################################################
# Assorted architecture dependent messages
###############################################################################

verb.unknown:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "PROBLEM: Dont know how to do this."
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexhp7:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
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
	@ echo "%M% %I% %E%"
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
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "not optimised sgi 64 bit compile for CC"
	@ echo ""
	@ echo "Note that this option is only for Matlab versions below 6.0"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.dll:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "Windows compile with gcc/mingw"
	@ echo "see http://www.mrc-cbu.cam.ac.uk/Imaging/gnumex20.html"
	@ echo "for instructions about installing gcc/mingw for"
	@ echo "compiling Mex files."
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexsol:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "Unix compile"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexlx:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "Linux compilation (Matlab 5.x) - using gcc"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexglx:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "Linux compilation (Matlab 6.x) - using gcc"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexa64:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "Linux compilation (Matlab 7.x, 64bits) - using gcc"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexaxp:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "keep your fingers crossed"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmac:
	@ echo "_____________________________________________________________"
	@ echo "%M% %I% %E%"
	@ echo ""
	@ echo "Unix compile for MacOS X"
	@ echo "_____________________________________________________________"
	@ echo ""

