@echo off
rem Windows compile with gcc using mingw
rem see http://www.mrc-cbu.cam.ac.uk/Imaging/gnumex20.html
rem for instructions about installing gcc for compiling Mex files.
rem Use the mingw gcc, and mexopts setting for mingw gcc, 
rem as described in the link above
rem 
rem Set mingw and matlab binary directories correctly here
set PATH=c:\mingw\bin;c:\matlabr11\bin;%PATH%
rem
rem Compiler directives
set DEFF=-DSPM_WIN32
set CC=gcc %DEFF%
set CMEX5=call mex.bat %DEFF%
rem 
echo SPM mex file compile for Windows
rem
echo Compiling SPM Windows utility files...
%CC% -c -o win32mmap.o win32mmap.c
%CMEX5% spm_win32utils.c
set added_objs=win32mmap.o spm_mapping.obj

echo %W% Matthew Brett %E%

echo Compiling volume utilities...
%CC% -c -o utils_uchar.o	spm_vol_utils.c -DSPM_UNSIGNED_CHAR
%CC% -c -o utils_short.o	spm_vol_utils.c -DSPM_SIGNED_SHORT
%CC% -c -o utils_int.o		spm_vol_utils.c -DSPM_SIGNED_INT
%CC% -c -o utils_schar.o	spm_vol_utils.c -DSPM_SIGNED_CHAR
%CC% -c -o utils_ushort.o	spm_vol_utils.c -DSPM_UNSIGNED_SHORT
%CC% -c -o utils_uint.o		spm_vol_utils.c -DSPM_UNSIGNED_INT
%CC% -c -o utils_float.o	spm_vol_utils.c -DSPM_FLOAT
%CC% -c -o utils_double.o	spm_vol_utils.c -DSPM_DOUBLE

rem Byteswapped images
%CC% -c -o utils_short_s.o	spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP
%CC% -c -o utils_int_s.o	spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP
%CC% -c -o utils_ushort_s.o	spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP
%CC% -c -o utils_uint_s.o	spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP
%CC% -c -o utils_float_s.o	spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP
%CC% -c -o utils_double_s.o	spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP

%CC% -c -o spm_make_lookup.o spm_make_lookup.c
%CC% -c -o spm_getdata.o spm_getdata.c
%CC% -c -o spm_vol_access.o  spm_vol_access.c 

rem utility routine
%CMEX5% -c spm_mapping.c

set vol_utils=utils_uchar.o utils_short.o utils_int.o utils_float.o utils_double.o utils_schar.o utils_ushort.o utils_uint.o utils_short_s.o utils_int_s.o utils_float_s.o utils_double_s.o utils_ushort_s.o utils_uint_s.o 	spm_vol_access.o spm_make_lookup.o spm_getdata.o %added_objs%

echo Adding to archive library spm_vol_utils.a...
del spm_vol_utils.a
ar rcv spm_vol_utils.a %vol_utils%
del %vol_utils%

echo Compiling mex files...
%CMEX5% spm_sample_vol.c	spm_vol_utils.a 
%CMEX5% spm_slice_vol.c		spm_vol_utils.a 
%CMEX5% spm_brainwarp.c		spm_vol_utils.a spm_matfuns.c
%CMEX5% spm_add.c		spm_vol_utils.a 
%CMEX5% spm_conv_vol.c		spm_vol_utils.a 
%CMEX5% spm_render_vol.c	spm_vol_utils.a 
%CMEX5% spm_global.c		spm_vol_utils.a 
%CMEX5% spm_resels_vol.c    	spm_vol_utils.a 
%CMEX5% spm_getxyz.c		spm_vol_utils.a 

%CMEX5% spm_atranspa.c
%CMEX5% spm_list_files.c
%CMEX5% spm_unlink.c
%CMEX5% spm_kronutil.c
%CMEX5% spm_project.c
%CMEX5% spm_hist2.c

%CMEX5% spm_max.c
%CMEX5% spm_clusters.c

echo Done.

