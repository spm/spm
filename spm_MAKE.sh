#!/bin/csh -vf
#
#     %W% John Ashburner %E%
#
# spm_MAKE will compile spm*.c scripts in a platform specific fashion
# see cmex
set CC = cc
set COPTIMFLAGS = "-xO5"
set cmex5 = "cmex5 COPTIMFLAGS=$COPTIMFLAGS"
set cmex4 = "cmex4 -O"

$CC -c -o utils_uchar.o spm_vol_utils.c -DUNSIGNED_CHAR $COPTIMFLAGS
$CC -c -o utils_short.o spm_vol_utils.c -DSIGNED_SHORT $COPTIMFLAGS
$CC -c -o utils_int.o spm_vol_utils.c -DSIGNED_INT $COPTIMFLAGS
$CC -c -o utils_float.o spm_vol_utils.c -DFLOAT $COPTIMFLAGS
$CC -c -o utils_double.o spm_vol_utils.c -DDOUBLE $COPTIMFLAGS

# Byteswapped images
$CC -c -o utils_short_s.o spm_vol_utils.c -DSIGNED_SHORT -DBYTESWAP $COPTIMFLAGS
$CC -c -o utils_int_s.o spm_vol_utils.c -DSIGNED_INT -DBYTESWAP $COPTIMFLAGS
$CC -c -o utils_float_s.o spm_vol_utils.c -DFLOAT -DBYTESWAP $COPTIMFLAGS
$CC -c -o utils_double_s.o spm_vol_utils.c -DDOUBLE -DBYTESWAP $COPTIMFLAGS

$CC -c -o spm_make_lookup.o spm_make_lookup.c -DDOUBLE -DBYTESWAP $COPTIMFLAGS
$CC -c -o spm_vol_utils2.o spm_vol_utils2.c $COPTIMFLAGS

ar rcv 	spm_vol_utils.a utils_uchar.o utils_short.o utils_int.o utils_float.o utils_double.o \
	utils_short_s.o utils_int_s.o utils_float_s.o utils_double_s.o \
	spm_vol_utils2.o spm_make_lookup.o
\rm utils_uchar.o utils_short.o utils_int.o utils_float.o utils_double.o \
	utils_short_s.o utils_int_s.o utils_float_s.o utils_double_s.o \
	spm_vol_utils2.o spm_make_lookup.o


$cmex5 spm_sample_vol.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_slice_vol.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_brainwarp.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_map_vol.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_unmap_vol.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_mean.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_conv_vol.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_render_vol.c	spm_vol_utils.a spm_mapping.c
$cmex5 spm_global.c	spm_vol_utils.a spm_mapping.c

$cmex5 spm_map_vol.c    spm_mapping.c
$cmex5 spm_unmap_vol.c  spm_mapping.c

$cmex5 spm_atranspa.c
$cmex4 spm_project.c
$cmex4 spm_max.c
$cmex4 spm_clusters.c
$cmex5 spm_list_files.c
$cmex5 spm_unlink.c
