#!make -f
#
# @(#)Makefile	1.3 John Ashburner 03/05/12
# $Id$
#
###############################################################################
#
# Suggestions for how to make this file a bit more elegant are welcome.  So far
# it works under SunOS and Linux (at the FIL) and some Windows systems
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

SunOS:
	make all SUF=mexsol  CC="cc -xO5"                   MEX="mex COPTIMFLAGS=-xO5"
Linux:
	make all SUF=mexglx  CC="gcc -O3 -funroll-loops"    MEX="mex COPTIMFLAGS='-O3 -funroll-loops'"
Linux.5:
	make all SUF=mexlx   CC="gcc -O3 -funroll-loops"    MEX="mex COPTIMFLAGS='-O3 -funroll-loops'"
HP-UX:
	make all SUF=mexhp7  CC="cc  -O +z -Ae +DAportable" MEX="mex COPTIMFLAGS=-O"
IRIX:
	make all SUF=mexsg   CC="cc  -O -mips2"             MEX="mex -O"
IRIX64:
	make all SUF=mexsg64 CC="cc  -O -mips4 -64"         MEX="mex -O"
AIX:
	make all SUF=mexrs6
OSF1:
	make all SUF=mexaxp
MAC:
	make all SUF=mexmac RANLIB="ranlib spm_vol_utils.mexmac.a"
#windows:
#	make all  SUF=dll MOSUF=obj CHMODIT="echo > null"  ADDED_OBS=win32mmap.dll.o ADDED_MEX=spm_win32utils.dll
windows:
	make all SUF=dll     CC="gcc -mno-cygwin -O3 -funroll-loops -DSPM_WIN32" MEX="mex.bat -DSPM_WIN32" MOSUF=obj CHMODIT="echo > NUL"

###############################################################################
# Architecture specific cleaning
###############################################################################
clean.SunOS:
	make clean SUF=mexsol
clean.Linux:
	make clean SUF=mexglx
clean.Linux.5:
	make clean SUF=mexlx
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

SPMMEX = spm_affdef.$(SUF) spm_def2det.$(SUF) spm_invdef.$(SUF)\
	 spm_invdef_vox.$(SUF) spm_warp.$(SUF)

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
	rm -f $(SPMMEX)

clean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting object (.o) files"
	@ echo "_____________________________________________________________"
	@ echo ""
	rm -f $(OBS)

###############################################################################
# Compile the mex files themselves
###############################################################################
MEX = mex -O
%.$(SUF) : %.c
	$(MEX) $<
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

verb.mexlx:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "Linux compilation (Matlab 5.x) - using gcc"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexaxp:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "keep your fingers crossed"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmac:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "Unix compile for MacOS X"
	@ echo "_____________________________________________________________"
	@ echo ""
