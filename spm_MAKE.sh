#!/bin/csh -v
#
#     %W% John Ashburner %E%
#
# spm_MAKE will compile spm*.c scripts in a platform specific fashion
# see cmex
alias CC    "cc -xO5"
alias cmex5 "cmex     COPTIMFLAGS=-xO5"
alias cmex4 "cmex -V4 COPTIMFLAGS=-xO5"

CC -c -o utils_uchar.o		spm_vol_utils.c -DUNSIGNED_CHAR
CC -c -o utils_short.o		spm_vol_utils.c -DSIGNED_SHORT
CC -c -o utils_int.o		spm_vol_utils.c -DSIGNED_INT
CC -c -o utils_schar.o		spm_vol_utils.c -DSIGNED_CHAR
CC -c -o utils_ushort.o	spm_vol_utils.c -DUNSIGNED_SHORT
CC -c -o utils_uint.o		spm_vol_utils.c -DUNSIGNED_INT
CC -c -o utils_float.o		spm_vol_utils.c -DFLOAT
CC -c -o utils_double.o	spm_vol_utils.c -DDOUBLE

# Byteswapped images
CC -c -o utils_short_s.o	spm_vol_utils.c -DSIGNED_SHORT -DBYTESWAP
CC -c -o utils_int_s.o		spm_vol_utils.c -DSIGNED_INT -DBYTESWAP
CC -c -o utils_ushort_s.o	spm_vol_utils.c -DUNSIGNED_SHORT -DBYTESWAP
CC -c -o utils_uint_s.o	spm_vol_utils.c -DUNSIGNED_INT -DBYTESWAP
CC -c -o utils_float_s.o	spm_vol_utils.c -DFLOAT -DBYTESWAP
CC -c -o utils_double_s.o	spm_vol_utils.c -DDOUBLE -DBYTESWAP

CC -c -o spm_make_lookup.o spm_make_lookup.c
CC -c -o spm_vol_utils2.o  spm_vol_utils2.c 

set vol_utils = (utils_uchar.o utils_short.o utils_int.o utils_float.o utils_double.o \
	utils_schar.o utils_ushort.o utils_uint.o \
	utils_short_s.o utils_int_s.o utils_float_s.o utils_double_s.o \
	utils_ushort_s.o utils_uint_s.o \
	spm_vol_utils2.o spm_make_lookup.o)

ar rcv spm_vol_utils.a $vol_utils
\rm $vol_utils

cmex5 spm_sample_vol.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_slice_vol.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_brainwarp.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_map_vol.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_unmap_vol.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_add.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_conv_vol.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_render_vol.c	spm_vol_utils.a spm_mapping.c
cmex5 spm_global.c	spm_vol_utils.a spm_mapping.c

cmex5 spm_atranspa.c
cmex5 spm_list_files.c
cmex5 spm_unlink.c

cmex5 spm_map_vol.c    spm_mapping.c
cmex5 spm_unmap_vol.c  spm_mapping.c
cmex4 spm_project.c
cmex4 spm_max.c
cmex4 spm_clusters.c
