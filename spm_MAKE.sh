#!/bin/sh
#
#     %W% JA, mods by Matthew Brett %E%
#
# spm_MAKE.sh will compile spm*.c scripts in a platform specific fashion
# see mex
#
# default compilation is for Unix
# pass string "windows" to script for Windows compile
# ie ./spm_MAKE.sh windows

if [ $# = 0 ]; then
	arch="unix";
else
	arch=$1;
fi

echo "SPM mex file compile for $arch"
echo 

case $arch in

    unix)
	# (default) unix compile
	#CC="gcc -Wall -O2"
	#cmex5="mex     COPTIMFLAGS=-O2"
	#cmex4="mex -V4 COPTIMFLAGS=-O2"
	CC="cc -xO5"
	cmex5="mex     COPTIMFLAGS=-xO5"
	cmex4="mex -V4 COPTIMFLAGS=-xO5"
	added_objs="spm_mapping.o";;
    windows)
	# set options for windows compile with gcc/mingw32
	#
	# see http://www.physiol.ox.ac.uk/~mb3/gnumex20.html
	# for instructions about installing gcc for
	# compiling Mex files.
	deff=-DSPM_WIN32
	CC="gcc -mno-cygwin -Wall $deff"
	cmex5="cmd /c mex $deff "
	cmex4="cmd /c mex $deff -V4 "

	# Windows added utility files
	$CC -c -o win32mmap.o win32mmap.c

	added_objs="win32mmap.o spm_mapping.obj";;

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

$cmex4 spm_max.c
$cmex4 spm_clusters.c

echo "Done."
