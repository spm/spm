#!/usr/bin/env make -f
# General Makefile to compile SPM C-MEX files
#
# Copyright (C) 1991-2023 Wellcome Centre for Human Neuroimaging
#
###############################################################################
#
# This Makefile has been tested under Linux, Windows and macOS.
#
# If you have to tweak this Makefile or Makefile.var to compile the SPM
# mex-files for your platform, please send the details to <fil.spm@ucl.ac.uk>
# so they can be included here.
#
# To compile and install SPM, type the following in a Terminal:
# >  make distclean && make && make install
#
# You can specify a particular platform with the following syntax:
# >  make PLATFORM=Your_Platform
# The standard targets are 'all', 'clean', 'distclean', 'doc', 'install',
# 'check' and 'test'.
#
# For a list of compatible compilers, see
#    https://www.mathworks.com/support/compilers.html
#
###############################################################################

include Makefile.var

###############################################################################
# Objects to go in the archive and mexfiles
###############################################################################

OBS     =\
	utils_uchar.$(MEXEXT).o utils_short.$(MEXEXT).o utils_int.$(MEXEXT).o \
	utils_schar.$(MEXEXT).o utils_ushort.$(MEXEXT).o utils_uint.$(MEXEXT).o\
	utils_float.$(MEXEXT).o utils_double.$(MEXEXT).o\
	utils_short_s.$(MEXEXT).o utils_int_s.$(MEXEXT).o\
	utils_ushort_s.$(MEXEXT).o utils_uint_s.$(MEXEXT).o\
	utils_float_s.$(MEXEXT).o utils_double_s.$(MEXEXT).o\
        utils_uint64.$(MEXEXT).o utils_int64.$(MEXEXT).o\
        utils_uint64_s.$(MEXEXT).o utils_int64_s.$(MEXEXT).o\
	spm_make_lookup.$(MEXEXT).o spm_getdata.$(MEXEXT).o spm_vol_access.$(MEXEXT).o\
	spm_mapping.$(MEXEXT).o

SPMMEX  =\
	spm_sample_vol.$(MEXEXT) spm_slice_vol.$(MEXEXT) spm_brainwarp.$(MEXEXT)\
	spm_conv_vol.$(MEXEXT) spm_render_vol.$(MEXEXT)\
	spm_global.$(MEXEXT) spm_resels_vol.$(MEXEXT)\
	spm_bsplinc.$(MEXEXT) spm_bsplins.$(MEXEXT)\
	spm_unlink.$(MEXEXT) spm_existfile.$(MEXEXT) spm_gamrnd.$(MEXEXT)\
	spm_hist.$(MEXEXT) spm_krutil.$(MEXEXT) spm_project.$(MEXEXT) spm_hist2.$(MEXEXT)\
	spm_dilate_erode.$(MEXEXT) spm_bwlabel.$(MEXEXT) spm_get_lm.$(MEXEXT)\
	spm_voronoi.$(MEXEXT) spm_mesh_dist.$(MEXEXT) spm_mesh_utils.$(MEXEXT) \
	spm_mrf.$(MEXEXT) spm_diffeo.$(MEXEXT) spm_field.$(MEXEXT) \
	spm_cat.$(MEXEXT) spm_jsonread.$(MEXEXT) spm_mesh_reduce.$(MEXEXT) \
	spm_mesh_geodesic.$(MEXEXT) spm_mesh_ray_triangle.${MEXEXT} spm_gmmlib.$(MEXEXT)

SUBDIRS =\
	@file_array/private \
	@gifti/private \
	@xmltree/private \
	toolbox/FieldMap

###############################################################################
# Public make targets
###############################################################################

all: verb.$(MEXEXT) main-all tbx-all verb.all.end

clean: verb.clean main-clean tbx-clean

distclean: clean verb.distclean main-distclean tbx-distclean verb.distclean.end

install: verb.install main-all tbx-all main-install tbx-install verb.install.end

external external-all external-clean external-distclean external-install external-tarball: verb.external
	@ $(MAKE) -C ../external --no-print-directory $(subst -,,$(subst external,,$@))

tarball: all tbx-tarball verb.tarball
	$(TAR) cf spm_mex.tar $(SPMMEX)
	@@for p in $(SUBDIRS); do \
		$(TAR) Af spm_mex.tar  ../$$p/spm_mex.tar;\
	done
	$(ZIP) spm_mex.tar

data:
	@ $(SPMEXE) eval "cd('../toolbox/mb');fil_install;"

check:
	@ $(SPMEXE) eval "spm_check_installation;"

test:
	@ $(SPMEXE) eval "spm_tests;close all force;"

###############################################################################
# Private make targets
###############################################################################

main-all: $(SPMMEX)

main-clean:
	$(DEL) $(OBS)

main-distclean:
	$(DEL) $(SPMMEX) spm_vol_utils.$(MEXEXT).a

main-install:
	$(COPY) $(SPMMEX) ..
	$(MOVE) ../spm_brainwarp.$(MEXEXT) ../toolbox/OldNorm

archive: spm_vol_utils.$(MEXEXT).a

tbx tbx-all tbx-clean tbx-distclean tbx-install tbx-tarball:
	@@for p in $(SUBDIRS); do \
		$(MAKE) -C ../$$p $(subst -,,$(subst tbx,,$@));\
	done

.PHONY: toolbox toolbox-all toolbox-clean toolbox-distclean toolbox-install toolbox-tarball
toolbox toolbox-all toolbox-clean toolbox-distclean toolbox-install toolbox-tarball:
	$(warning Target "toolbox" is deprecated.)

.NOTPARALLEL:

###############################################################################
# Compile spm_vol_utils.c with various flags
###############################################################################

spm_vol_utils.$(MEXEXT).a: $(OBS)
	$(DEL) $@
ifdef ARNOSPACE
	$(AR)$@ $(OBS)
else
	$(AR) $@ $(OBS)
endif

UTILS = spm_vol_utils.c spm_make_lookup.h spm_getdata.h

utils_uchar.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_CHAR $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_short.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_SHORT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_int.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_INT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_schar.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_CHAR $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_ushort.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_uint.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_INT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_float.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_FLOAT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_double.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_DOUBLE $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_uint64.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_LONG_LONG $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_int64.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_LONG_LONG $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_short_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_int_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_ushort_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_uint_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_float_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_double_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_uint64_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_LONG_LONG -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_int64_s.$(MEXEXT).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_LONG_LONG -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

###############################################################################
# Compile a few additional C routines for linking
###############################################################################

%.$(MEXEXT).o : %.c
	$(MEX) -c $< $(MEXEND)
	$(MOVE) %.$(MOSUF) $@

spm_getdata.$(MEXEXT).o: spm_getdata.c spm_getdata.h
	$(MEX) -c spm_getdata.c $(MEXEND)
	$(MOVE) spm_getdata.$(MOSUF) $@
	
spm_vol_access.$(MEXEXT).o: spm_vol_access.c spm_vol_access.h spm_datatypes.h
	$(MEX) -c spm_vol_access.c $(MEXEND)
	$(MOVE) spm_vol_access.$(MOSUF) $@

spm_make_lookup.$(MEXEXT).o: spm_make_lookup.c spm_make_lookup.h
	$(MEX) -c spm_make_lookup.c $(MEXEND)
	$(MOVE) spm_make_lookup.$(MOSUF) $@
	
spm_mapping.$(MEXEXT).o: spm_mapping.c spm_mapping.h spm_vol_access.h spm_datatypes.h
	$(MEX) -c spm_mapping.c $(MEXEND)
	$(MOVE) spm_mapping.$(MOSUF) $@

###############################################################################
# Compile the mex files themselves
###############################################################################

%.$(MEXEXT) : %.c
	$(MEX) $< $(MEXEND)

spm_brainwarp.$(MEXEXT): spm_brainwarp.c spm_matfuns.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_brainwarp.c spm_vol_utils.$(MEXEXT).a spm_matfuns.c $(MEXEND)

spm_bsplinc.$(MEXEXT): spm_bsplinc.c bsplines.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h bsplines.h
	$(MEX) spm_bsplinc.c bsplines.c spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_conv_vol.$(MEXEXT): spm_conv_vol.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h spm_datatypes.h
	$(MEX) spm_conv_vol.c spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_global.$(MEXEXT): spm_global.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_global.c spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_render_vol.$(MEXEXT): spm_render_vol.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_render_vol.c spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_resels_vol.$(MEXEXT): spm_resels_vol.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_resels_vol.c spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_sample_vol.$(MEXEXT): spm_sample_vol.c spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_sample_vol.c spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_slice_vol.$(MEXEXT): spm_slice_vol.c  spm_vol_utils.$(MEXEXT).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_slice_vol.c  spm_vol_utils.$(MEXEXT).a $(MEXEND)

spm_bsplins.$(MEXEXT): spm_bsplins.c bsplines.c bsplines.h
	$(MEX) spm_bsplins.c bsplines.c $(MEXEND)
	
spm_bwlabel.$(MEXEXT): spm_bwlabel.c
	$(MEX) spm_bwlabel.c $(MEXEND)

spm_dilate_erode.$(MEXEXT): spm_dilate_erode.c
	$(MEX) spm_dilate_erode.c $(MEXEND)

spm_existfile.$(MEXEXT): spm_existfile.c
	$(MEX) spm_existfile.c $(MEXEND)

spm_gamrnd.$(MEXEXT): spm_gamrnd.c
	$(MEX) spm_gamrnd.c $(MEXEND)

spm_get_lm.$(MEXEXT): spm_get_lm.c
	$(MEX) spm_get_lm.c $(MEXEND)

spm_hist.$(MEXEXT): spm_hist.c
	$(MEX) spm_hist.c $(MEXEND)

spm_hist2.$(MEXEXT): spm_hist2.c hist2.c hist2.h
	$(MEX) spm_hist2.c hist2.c $(MEXEND)

spm_krutil.$(MEXEXT): spm_krutil.c
	$(MEX) spm_krutil.c $(MEXEND)

spm_project.$(MEXEXT): spm_project.c
	$(MEX) spm_project.c $(MEXEND)

spm_unlink.$(MEXEXT): spm_unlink.c
	$(MEX) spm_unlink.c $(MEXEND)

spm_voronoi.$(MEXEXT): spm_voronoi.c
	$(MEX) spm_voronoi.c $(MEXEND)

spm_mrf.$(MEXEXT): spm_mrf.c
	$(MEX) spm_mrf.c $(MEXEND)

spm_diffeo.$(MEXEXT): spm_diffeo.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c\
		shoot_regularisers.c shoot_expm3.c shoot_invdef.c shoot_dartel.c\
		shoot_boundary.c spm_openmp.c shoot_bsplines.c bsplines.c
	$(MEX) spm_diffeo.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c\
		shoot_regularisers.c shoot_expm3.c shoot_invdef.c shoot_dartel.c\
		shoot_boundary.c spm_openmp.c shoot_bsplines.c bsplines.c\
		-DIMAGE_SINGLE $(MEXEND)

spm_field.$(MEXEXT): spm_field.c  shoot_optimN.c shoot_multiscale.c shoot_boundary.c\
		spm_openmp.c
	$(MEX) spm_field.c shoot_optimN.c shoot_multiscale.c shoot_boundary.c\
		spm_openmp.c $(MEXEND)

spm_gmmlib.$(MEXEXT): spm_gmmlib.c gmmlib.c gmmlib.h 
	$(MEX) spm_gmmlib.c gmmlib.c -largeArrayDims $(MEXEND)

spm_jsonread.$(MEXEXT): spm_jsonread.c external/jsmn/jsmn.h external/jsmn/jsmn.c
	$(MEX) spm_jsonread.c external/jsmn/jsmn.c -DJSMN_PARENT_LINKS $(MEXEND)

spm_mesh_dist.$(MEXEXT):spm_mesh_dist.cpp external/TriangleMeshDistance/TriangleMeshDistance.h
	$(MEX) spm_mesh_dist.cpp $(MEXEND)

spm_mesh_geodesic.$(MEXEXT): spm_mesh_geodesic.cpp
	$(MEX) spm_mesh_geodesic.cpp $(MEXEND)

spm_mesh_reduce.$(MEXEXT): spm_mesh_reduce.c external/nii2mesh/quadric.c\
		external/nii2mesh/quadric.h external/nii2mesh/meshtypes.h
	$(MEX) spm_mesh_reduce.c external/nii2mesh/quadric.c $(MEXEND)

spm_mesh_ray_triangle.$(MEXEXT): spm_mesh_ray_triangle.c external/raytri/raytri.c
	$(MEX) spm_mesh_ray_triangle.c $(MEXEND)


###############################################################################
# Compile shared library libSPM
###############################################################################

libSPM.so: shoot_diffeo3d.c shoot_optimN.c shoot_boundary.c shoot_multiscale.c\
	   shoot_regularisers.c shoot_optim3d.c hist2.c shoot_invdef.c shoot_expm3.c\
	   gmmlib.c spm_openmp.c\
	   shoot_diffeo3d.h shoot_optimN.h shoot_boundary.h shoot_multiscale.h\
	   shoot_regularisers.h shoot_optim3d.h hist2.h shoot_expm3.h gmmlib.h\
	   spm_openmp.h
	$(CC) -O4 -shared -fPIC -Wl,-soname,$@ -o $@ -lm -DIMAGE_SINGLE\
	   shoot_diffeo3d.c shoot_optimN.c shoot_boundary.c shoot_multiscale.c\
	   shoot_regularisers.c shoot_optim3d.c hist2.c shoot_invdef.c shoot_expm3.c\
	   gmmlib.c spm_openmp.c  $(OMPFLAG)

###############################################################################
# Display Messages
###############################################################################

verb.clean:
	$(call verb, "Deleting object (.o) files")

verb.distclean:
	$(call verb, "Deleting MEX (.$(MEXEXT)) and archive (.a) files")

verb.install:
	$(call verb, "Installing MEX files")

verb.tarball:
	$(call verb, "Creating archive spm_mex.tar.gz")

verb.mexw64:
	$(call verb, "Windows compilation (64 bit)")

verb.mexa64:
	$(call verb, "Linux compilation (x86-64)")

verb.mexmaci64:
	$(call verb, "macOS compilation (Intel 64 bit)")

verb.mexmaca64:
	$(call verb, "macOS compilation (ARM 64 bit)")

verb.mex:
	$(call verb, "${PLATFORM} compilation ($(shell ${MEXBIN} -v | head -n 1))")

verb.all.end:
	$(call verb, "Compilation: done")

verb.distclean.end:
	$(call verb, "Distclean: done")

verb.install.end:
	$(call verb, "Installation: done")

verb.external:
	$(call verb, "In external")
