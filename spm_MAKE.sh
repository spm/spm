#!/bin/sh
#
#     %W% JA, mods by Matthew Brett %E%
#
# spm_MAKE.sh will compile spm*.c scripts in a platform specific fashion
# see mex
#
# default compilation is for Sun cc
# pass string to script from list below for other compile types
# e.g. ./spm_MAKE.sh windows

if [ $# = 0 ]; then
	arch="sun";
else
	if [ $1 = "--help" ]; then
		echo "spm_MAKE [architecture/compiler]"
		echo "Call with architecture/compiler as argument, where"
		echo "architecture/compiler may be:"
		echo "   sun (default)   - Sun (SunOS, solaris) ?DEC, using cc"
		echo "   windows         - windows (NT, 95/98), using EGCS gcc"
		echo "   gcc             - gcc compile for unix, including Sun, linux"
		echo "   sgi             - Irix 32 bit compile with cc"
		echo "   sgi64           - Irix 64 bit compile with cc"
		echo "   hpux            - ?HPUX (HP) or ??AIX (IBM) compile with cc"
		exit
	else
		arch=$1;
	fi
fi

echo
echo "SPM mex file compile for $arch"
echo 

added_objs="spm_mapping.o"

case $arch in
    sun)
	# (default) unix compile for Sun CC
	CC="cc -xO5"
	cmex5="mex     COPTIMFLAGS=-xO5"
	cmex4="mex -V4 COPTIMFLAGS=-xO5";;
    windows)
	# windows compile with EGCS gcc/mingw32
	# see http://www.physiol.ox.ac.uk/~mb3/gnumex20.html
	# for instructions about installing gcc for
	# compiling Mex files.
	deff=-DSPM_WIN32
	CC="gcc -mno-cygwin $deff"
	cmex5="mex.bat $deff "
	cmex4="mex.bat $deff -V4 "
	# Windows added utility files
	$CC -c -o win32mmap.o win32mmap.c
	$cmex5 spm_win32utils.c
	added_objs="win32mmap.o spm_mapping.obj";;
    gcc)
	# optimised standard unix compile for gcc
	# this should work on Sun, Linux etc
	CC="gcc -O2"
	cmex5="mex     COPTIMFLAGS=-O2"
	cmex4="mex -V4 COPTIMFLAGS=-O2";;
    sgi)
	# not optimised unix compile for CC
	CC="cc"
	cmex5="mex"
	cmex4="mex -V4";;
    sgi64)
	# not optimised sgi 64 bit compile for CC
	CC="cc -64"
	cmex5="mex"
	cmex4="mex -V4";;
    hpux)
	# unix compile for hpux cc, and maybe aix cc
	CC="cc -O +Z"
	cmex5="mex     COPTIMFLAGS=-O"
	cmex4="mex -V4 COPTIMFLAGS=-O";;
   *)
	echo "Sorry, not set up for architecture $arch"
	exit;;
esac

echo "Compiling volume utilities..."
$CC -c -o utils_uchar.o		spm_vol_utils.c -DSPM_UNSIGNED_CHAR
$CC -c -o utils_short.o		spm_vol_utils.c -DSPM_SIGNED_SHORT
$CC -c -o utils_int.o		spm_vol_utils.c -DSPM_SIGNED_INT
$CC -c -o utils_schar.o		spm_vol_utils.c -DSPM_SIGNED_CHAR
$CC -c -o utils_ushort.o	spm_vol_utils.c -DSPM_UNSIGNED_SHORT
$CC -c -o utils_uint.o		spm_vol_utils.c -DSPM_UNSIGNED_INT
$CC -c -o utils_float.o		spm_vol_utils.c -DSPM_FLOAT
$CC -c -o utils_double.o	spm_vol_utils.c -DSPM_DOUBLE

# Byteswapped images
$CC -c -o utils_short_s.o	spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP
$CC -c -o utils_int_s.o		spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP
$CC -c -o utils_ushort_s.o	spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP
$CC -c -o utils_uint_s.o	spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP
$CC -c -o utils_float_s.o	spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP
$CC -c -o utils_double_s.o	spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP

$CC -c -o spm_make_lookup.o spm_make_lookup.c
$CC -c -o spm_getdata.o spm_getdata.c
$CC -c -o spm_vol_access.o  spm_vol_access.c 

# utility routine
$cmex5 -c spm_mapping.c

vol_utils="utils_uchar.o utils_short.o utils_int.o utils_float.o utils_double.o \
	utils_schar.o utils_ushort.o utils_uint.o \
	utils_short_s.o utils_int_s.o utils_float_s.o utils_double_s.o \
	utils_ushort_s.o utils_uint_s.o \
	spm_vol_access.o spm_make_lookup.o spm_getdata.o $added_objs"

echo "Adding to archive library spm_vol_utils.a..."
\rm spm_vol_utils.a
ar rcv spm_vol_utils.a $vol_utils
\rm $vol_utils

echo "Compiling mex files..."
$cmex5 spm_sample_vol.c	spm_vol_utils.a 
$cmex5 spm_slice_vol.c	spm_vol_utils.a 
$cmex5 spm_brainwarp.c	spm_vol_utils.a 
$cmex5 spm_map_vol.c	spm_vol_utils.a 
$cmex5 spm_unmap_vol.c	spm_vol_utils.a 
$cmex5 spm_add.c	spm_vol_utils.a 
$cmex5 spm_conv_vol.c	spm_vol_utils.a 
$cmex5 spm_render_vol.c	spm_vol_utils.a 
$cmex5 spm_global.c	spm_vol_utils.a 
$cmex5 spm_resels_vol.c	spm_vol_utils.a 
$cmex5 spm_getxyz.c	spm_vol_utils.a 

$cmex5 spm_atranspa.c
$cmex5 spm_list_files.c
$cmex5 spm_unlink.c
$cmex5 spm_kronutil.c
$cmex5 spm_project.c

$cmex5 spm_max.c
$cmex5 spm_clusters.c

echo "Done."
