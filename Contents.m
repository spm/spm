% SPM96
% Statistical Parametric Mapping - SPM96
%_______________________________________________________________________
% 
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  University College London
%
% John Ashburner, Karl Friston, Andrew Holmes, Jean-Baptiste Poline
%_______________________________________________________________________
%
% SPM96 - Manifest
%
% SPM96 is written for Matlab 4.2c under UNIX
%-----------------------------------------------------------------------
% README                   1.1      (New for SPM96)
% Contents.m               %I%
% Grid.mat
% MIP.mat                  1.2
% Split.mat
% connex.h                 1.2
% dbh.h                    (Mayo clinic database sub-definitions)
% render.mat
% spm.m                    1.26     (1.25)
% spm.man                  1.19     (1.17)
% spm_AnCova.m             1.4
% spm_DesMtx.m             1.1
% spm_DesMtxSca.m          1.1
% spm_F.man                1.2
% spm_Fcdf.m               1.2
% spm_Fpdf.m               1.2
% spm_Gcdf.m               1.2
% spm_Gpdf.m               1.2
% spm_MAKE
% spm_Ncdf.m               1.3
% spm_Npdf.m               1.3
% spm_P.m                  1.2
% spm_Pcdf.m               1.1
% spm_Pkn.m                1.3
% spm_Pn.m                 1.1
% spm_Ppdf.m               1.1
% spm_Pz.m                 1.1
% spm_Tcdf.m               1.2
% spm_Tpdf.m               1.3
% spm_W.m                  1.2
% spm_Xcdf.m               1.2
% spm_Xpdf.m               1.2
% spm_affsub1.m            1.1
% spm_affsub2.m            1.4
% spm_affsub3.m            1.5
% spm_append.m             1.1
% spm_atranspa.c           1.1
% spm_atranspa.m           1.2
% spm_atranspa.mexsg
% spm_atranspa.mexsol
% spm_average.m            1.4
% spm_bb.m                 1.1
% spm_box.c                1.1
% spm_box.m                1.1
% spm_box.mexsg
% spm_box.mexsol
% spm_brainwarp.c          1.8
% spm_brainwarp.m          1.4
% spm_brainwarp.mexsg
% spm_brainwarp.mexsol
% spm_choose.m             1.1
% spm_clf.m                1.2
% spm_clusters.c           1.1
% spm_clusters.m           1.1
% spm_clusters.mexsg
% spm_clusters.mexsol
% spm_conv.m               1.3      (1.2 )
% spm_conv_vol.c           1.3
% spm_conv_vol.m           1.2
% spm_conv_vol.mexsg
% spm_conv_vol.mexsol
% spm_coregister.m         1.12     (1.11)
% spm_dctmtx.m             1.3
% spm_defaults.m           1.8      (1.7 )
% spm_defaults_edit.m      1.9
% spm_detrend.m            1.2
% spm_display.m            1.5
% spm_en.m                 1.1
% spm_figure.m             1.14     (1.13)
% spm_fix_header.m         1.1
% spm_fmri.man             1.2
% spm_fmri_spm_ui.m        1.17
% spm_format.man           1.2
% spm_fzero.m              1.1
% spm_get.m                1.13     (1.12)
% spm_get_space.m          1.7      (1.6 )
% spm_global.c             1.1
% spm_global.m             1.2
% spm_global.mexsg
% spm_global.mexsol
% pm_graph.m               1.3
% spm_grid.m               1.1
% spm_header_edit.m        1.2
% spm_help.m               1.15
% spm_hread.m              1.4
% spm_hwrite.m             1.4
% spm_image.m              1.4
% spm_image.man            1.1
% spm_image_funks.m        1.4
% spm_imatrix.m            1.2       (New for SPM96)
% spm_input.m              1.9
% spm_invFcdf.m            1.2
% spm_invGcdf.m            1.2
% spm_invNcdf.m            1.2
% spm_invTcdf.m            1.2
% spm_invXcdf.m            1.2
% spm_invkcdf.m            1.1
% spm_k.m                  1.1
% spm_kcdf.m               1.1
% spm_kpdf.m               1.1
% spm_lambda.m             1.6
% spm_lctx.mat
% spm_list_files.c         1.6
% spm_list_files.m         1.2
% spm_list_files.mexsg
% spm_list_files.mexsol
% spm_load.m               1.3
% spm_log.m                1.1
% spm_map.m                1.3
% spm_map_vol.c            1.7
% spm_map_vol.m            1.1
% spm_map_vol.mexsg
% spm_map_vol.mexsol
% spm_mat.man              1.1
% spm_matrix.m             1.1
% spm_max.c                1.1
% spm_max.m                1.1
% spm_max.mexsg
% spm_max.mexsol
% spm_maxima.m             1.8
% spm_mean.c               1.10
% spm_mean.m               1.3
% spm_mean.mexsg
% spm_mean.mexsol
% spm_meanby.m             1.5
% spm_methods.man          1.1
% spm_min_Pn.m             1.1
% spm_min_Pz.m             1.1
% spm_mip.m                1.8
% spm_mip_ui.m             1.3
% spm_modality.man         1.1
% spm_motd.man             1.2      (1.1 )
% spm_orthviews.m          1.2
% spm_pF.m                 1.2
% spm_pet.man              1.2
% spm_picture.m            1.2
% spm_print.m              1.3
% spm_progress_bar.m       1.4
% spm_project.c            1.6
% spm_project.m            1.1
% spm_project.mexsg
% spm_project.mexsol
% spm_projections.m        1.12     (1.11)
% spm_projections.man      1.1
% spm_projectionsF.m       1.3
% spm_projectionsF_ui.m    1.6
% spm_projections_ui.m     1.11.1.1 (1.11)
% spm_rctx.mat
% spm_readXA.m             1.1
% spm_realign.m            1.16     (1.15)
% spm_realign.man          1.3
% spm_render.m             1.4
% spm_render_vol.c         1.3
% spm_render_vol.m         1.1
% spm_render_vol.mexsg
% spm_render_vol.mexsol
% spm_renviews.m           1.1
% spm_resize.m             1.1
% spm_results.man          1.1
% spm_results_ui.m         1.7
% spm_sample_vol.c         1.8
% spm_sample_vol.m         1.1
% spm_sample_vol.mexsg
% spm_sample_vol.mexsol
% spm_sections.m           1.7
% spm_segment.m            1.27     (1.25)
% spm_slice_vol.c          1.5
% spm_slice_vol.m          1.1
% spm_slice_vol.mexsg
% spm_slice_vol.mexsol
% spm_smooth.m             1.4
% spm_smooth.man           1.1
% spm_smooth_ui.m          1.3
% spm_sn3d.m               1.16     (1.15) 
% spm_sn3d.man             1.2
% spm_snbasis_map.m        1.5
% spm_spm.m                1.23
% spm_spm.man              1.1
% spm_spm_ui.m             1.15.1.2 (1.15)
% spm_sptop.m              1.2      (1.1 )
% spm_str_manip.m          1.3
% spm_svd.m                1.3
% spm_svd.man              1.1
% spm_svd_ui.m             1.2
% spm_t2z.m                1.8
% spm_transverse.m         1.7
% spm_type.m               1.2      (1.1 )
% spm_unlink.c             1.1
% spm_unlink.m             1.2
% spm_unlink.mexsg
% spm_unlink.mexsol
% spm_unmap.m              1.1
% spm_unmap_vol.c          1.2
% spm_unmap_vol.m          1.1
% spm_unmap_vol.mexsg
% spm_unmap_vol.mexsol
% spm_vol_utils.c          1.5
% spm_write.m              1.2
% spm_write_filtered.m     1.1
% spm_write_sn.m           1.5      (1.4 )
% templates.man            1.1
% volume.h                 1.6
% 
% apriori
% apriori/csf.hdr
% apriori/csf.img
% apriori/gray.hdr
% apriori/gray.img
% apriori/white.hdr
% apriori/white.img
% apriori/symmetric_csf.hdr
% apriori/symmetric_csf.img
% apriori/symmetric_gray.hdr
% apriori/symmetric_gray.img
% apriori/symmetric_white.hdr
% apriori/symmetric_white.img
% 
% canonical
% canonical/avg305T1.hdr
% canonical/avg305T1.img
% canonical/T1.img
% canonical/T1.hdr
% 
% templates
% templates/PET.hdr
% templates/PET.img
% templates/T1.hdr
% templates/T1.img
% templates/T2.hdr
% templates/T2.img
%
%-----------------------------------------------------------------------
% The following SPM96 files have been updated since the beta release:
%
% README                1.1             - Basic README file (New for SPM96)
% Contents.m            %I%             - This file, Manifest updated
% spm_motd.man          1.2      (1.1 ) - Message of the day file updated
% spm.man               1.19     (1.17) - Truncated paragraph, FLIP explained
% spm_type.m            1.2      (1.1 ) + Spelling corrections!
% spm.m                 1.26     (1.25) ~ Watermark removal
% spm_figure.m          1.14     (1.13) ~ Bug fix: redundent findobj
% spm_get.m             1.13     (1.12) ~ Bug in CmdLine code fixed
% spm_spm_ui.m          1.15.1.2 (1.15) ~ Explicit constant, cCovNoInt int's
% spm_projections_ui.m  1.11.1.1 (1.11) ~ svd problems(ML4.2c/Sol2.4)
% spm_projections.m     1.12     (1.11) ~ Bug fix: Printed FWHMvoxels
% spm_sptop.m           1.2      (1.1 ) ~ Kernel normalisation
% spm_conv.m            1.3      (1.2 ) + Uses new spm_sptop.m
% spm_realign.m         1.16     (1.15) ~ Old (x,y&z) adjustment reinstated
% spm_imatrix.m         1.2             ~ New file for spm_realign.m v1.16
% spm_coregister.m      1.12     (1.11) ~ Deletion of temporary files
% spm_defaults.m        1.8      (1.7 ) + sptl_Rglrztn addded
% spm_sn3d.m            1.16     (1.15) + Global sptl_Rglrztn
% spm_get_space.m       1.7      (1.6 ) ~ VX : voxel sizes occasionally
% spm_segment.m         1.27     (1.25) + niter=48, default starting estimate
% spm_write_sn.m        1.5      (1.4 ) ~ Unmap files
%
% Version numbers in brackets refer to the version distributed with
% SPM96b. Files labelled "~" are publicised updates to SPM96b,
% discussed on the help list <spm@mailbase.ac.uk> and on the SPMweb
% site http://www.fil.ion.ucl.ac.uk/spm. Files labelled "+" are
% cosmetic improvements over SPM96b.
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%=======================================================================
% PROGRAMMERS NOTE:
% This (Contents.m) is the contents file for SPM, used by spm('Ver') to
% recover the version number. MatLab's ver also uses Contents.m files
% to identify toolbox versions.
% Line1: Version (one word). Line2: One line description
%=======================================================================
