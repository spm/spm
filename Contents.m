% SPM99b (c) 1991,1994-1999 : The Wellcome Department of Cognitive Neurology
% Statistical Parametric Mapping       -       SPM99 public beta release
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  SPM - http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________
%
% This Contents.m file holds the version ID for this release of Matlab,
% and contains a manifest of the included functions and their version numbers.
%
% SPM99 is written for Matlab v5.2.1 under UNIX and Windows
% ( Compiled binaries of external MEX functions are provided for:       )
% (                   Solaris2, Linux, and Windows                     )
%
% See spm.man for details of this release.
% See the README for information on installation and getting started.
% See spm_motd.man for last minute release details.
%
% SPM (being the collection of files given in the manifest below) is
% free but copyright software, distributed under the terms of the GNU
% General Public Licence as published by the Free Software Foundation
% (either version 2, as given in file spm_LICENCE.man, or at your option,
% any later version). Further details on "copyleft" can be found at
% http://www.gnu.org/copyleft/.
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%
%
% SPM99b - Manifest
%-----------------------------------------------------------------------
% Contents.m            	2.2
% Grid.mat                      1.1
% MIP.mat                       1.2
% README.txt                    2.4
% Split.mat                     1.1
% connex.h                      1.2
% dbh.h                         (Mayo clinic database sub-definitions)
% spm.m                         2.51
% spm.man                       2.5
% spm_Bcdf.m                    2.2
% spm_Bpdf.m                    2.2
% spm_DesMtx.m                  2.9
% spm_DesRep.m                  2.13
% spm_FcUtil.m                  2.8
% spm_Fcdf.m                    2.2
% spm_FcnTpl.txt                1.1
% spm_Fpdf.m                    2.2
% spm_Gcdf.m                    2.2
% spm_Gpdf.m                    2.2
% spm_Icdf.m                    2.2
% spm_Ipdf.m                    2.2
% spm_LICENCE.man               2.4
% spm_MAKE.sh                   2.14
% spm_Ncdf.m                    2.2
% spm_Npdf.m                    2.2
% spm_P.m                       2.2
% spm_Pcdf.m                    2.2
% spm_Pec_resels.m              2.1
% spm_Ppdf.m                    2.2
% spm_RandFX.man                1.4
% spm_SpUtil.m                  2.12
% spm_Tcdf.m                    2.2
% spm_Tpdf.m                    2.2
% spm_VOI.m                     2.8
% spm_Vintrinsic.m              2.2
% spm_Volterra.m                2.1
% spm_XYZreg.m                  2.4
% spm_XYZreg_Ex1.m              2.1
% spm_XYZreg_Ex2.m              2.1
% spm_Xcdf.m                    2.2
% spm_Xpdf.m                    2.2
% spm_add.c                     2.8
% spm_add.dll                   2.8
% spm_add.m                     2.7
% spm_add.mexlx                 2.8
% spm_add.mexsol                2.8
% spm_adjmean_fmri_ui.m         2.6
% spm_adjmean_ui.m              2.5
% spm_affsub3.m                 2.5
% spm_append.m                  2.2
% spm_atranspa.c                2.1
% spm_atranspa.dll              2.1
% spm_atranspa.m                2.1
% spm_atranspa.mexlx            2.1
% spm_atranspa.mexsol           2.1
% spm_bigend.m                  2.1
% spm_brainwarp.c               2.8
% spm_brainwarp.dll             2.8
% spm_brainwarp.m               2.2
% spm_brainwarp.mexlx           2.8
% spm_brainwarp.mexsol          2.8
% spm_check_registration.m      1.1
% spm_chi2_plot.m               1.2
% spm_choose.m                  2.6
% spm_clf.m                     2.1
% spm_clusters.c                2.2
% spm_clusters.dll              2.2
% spm_clusters.m                2.2
% spm_clusters.mexlx            2.2
% spm_clusters.mexsol           2.2
% spm_conman.m                  2.11
% spm_conv.m                    1.4
% spm_conv_vol.c                2.3
% spm_conv_vol.dll              2.3
% spm_conv_vol.m                2.1
% spm_conv_vol.mexlx            2.3
% spm_conv_vol.mexsol           2.3
% spm_coregister.m              2.5
% spm_create_image.m            2.4
% spm_datatypes.h               2.1
% spm_dbm.m                     2.2
% spm_dctmtx.m                  1.3
% spm_defaults.m                2.9
% spm_defaults_edit.m           2.5
% spm_detrend.m                 2.1
% spm_dummy.m                   2.1
% spm_eigenim.m                 2.1
% spm_en.m                      1.1
% spm_extract.m                 2.2
% spm_fMRI_design.m             2.16
% spm_fMRI_design_show.m        2.13
% spm_figure.m                  2.21
% spm_filter.m                  2.1
% spm_fir.m                     2.1
% spm_fmri.man                  1.4
% spm_fmri_spm_ui.m             2.25
% spm_format.man                2.1
% spm_get.m                     2.30
% spm_getSPM.m                  2.20
% spm_get_bf.m                  2.10
% spm_get_ons.m                 2.15
% spm_get_space.m               2.4
% spm_getdata.c                 2.3
% spm_getdata.h                 2.2
% spm_getxyz.c                  2.3
% spm_getxyz.dll                2.3
% spm_getxyz.m                  2.2
% spm_getxyz.mexlx              2.3
% spm_getxyz.mexsol             2.3
% spm_global.c                  2.4
% spm_global.dll                2.4
% spm_global.m                  2.2
% spm_global.mexlx              2.4
% spm_global.mexsol             2.4
% spm_graph.m                   2.19
% spm_grid.m                    1.1
% spm_header_edit.m             2.1
% spm_help.m                    2.32
% spm_hread.m                   2.3
% spm_hrf.m                     2.7
% spm_hwrite.m                  2.1
% spm_image.m                   2.14
% spm_imatrix.m                 2.1
% spm_imcalc.m                  2.6
% spm_imcalc_ui.m               2.6
% spm_input.m                   2.38
% spm_invBcdf.m                 2.2
% spm_invFcdf.m                 2.2
% spm_invGcdf.m                 2.2
% spm_invIcdf.m                 2.2
% spm_invNcdf.m                 2.2
% spm_invPcdf.m                 2.2
% spm_invTcdf.m                 2.2
% spm_invXcdf.m                 2.3
% spm_kronutil.c                2.1
% spm_kronutil.dll              2.1
% spm_kronutil.m                2.2
% spm_kronutil.mexlx            2.1
% spm_kronutil.mexsol           2.1
% spm_list.m                    2.21
% spm_list_files.c              2.7
% spm_list_files.dll            2.7
% spm_list_files.m              2.1
% spm_list_files.mexlx          2.7
% spm_list_files.mexsol         2.7
% spm_load.m                    2.1
% spm_log.m                     2.2
% spm_make_filter.m             2.4
% spm_make_lookup.c             2.3
% spm_make_lookup.h             2.1
% spm_map.h                     2.2
% spm_map.m                     2.1
% spm_map_vol.c                 2.1
% spm_map_vol.dll               2.1
% spm_map_vol.m                 2.1
% spm_map_vol.mexlx             2.1
% spm_map_vol.mexsol            2.1
% spm_mapping.c                 2.3
% spm_mapping.h                 2.2
% spm_mask.m                    2.7
% spm_matfuns.c                 2.1
% spm_matrix.m                  1.1
% spm_matx.m                    2.4
% spm_max.c                     2.3
% spm_max.dll                   2.3
% spm_max.m                     2.2
% spm_max.mexlx                 2.3
% spm_max.mexsol                2.3
% spm_maxima.m                  2.8
% spm_mean_ui.m                 2.4
% spm_meanby.m                  2.1
% spm_mip.m                     2.3
% spm_mip_ui.m                  2.8
% spm_motd.man                  2.2
% spm_mvNpdf.m                  2.1
% spm_nCr.m                     2.1
% spm_orthviews.m               2.16
% spm_pet.man                   1.4
% spm_platform.m                2.4
% spm_print.m                   1.3
% spm_progress_bar.m            2.1
% spm_project.c                 2.2
% spm_project.dll               2.2
% spm_project.m                 2.2
% spm_project.mexlx             2.2
% spm_project.mexsol            2.2
% spm_read_vols.m               2.4
% spm_realign.m                 2.16
% spm_render.m                  2.15
% spm_render_vol.c              2.2
% spm_render_vol.dll            2.2
% spm_render_vol.m              2.1
% spm_render_vol.mexlx          2.2
% spm_render_vol.mexsol         2.2
% spm_renviews.m                2.1
% spm_resels.m                  2.1
% spm_resels_vol.c              2.4
% spm_resels_vol.dll            2.4
% spm_resels_vol.m              2.2
% spm_resels_vol.mexlx          2.4
% spm_resels_vol.mexsol         2.4
% spm_resize.m                  1.1
% spm_resss.m                   2.7
% spm_results.m                 2.7
% spm_results_ui.m              2.22
% spm_sample_vol.c              2.1
% spm_sample_vol.dll            2.1
% spm_sample_vol.m              2.1
% spm_sample_vol.mexlx          2.1
% spm_sample_vol.mexsol         2.1
% spm_sections.m                2.12
% spm_segment.m                 2.9
% spm_slice_timing.m            2.2
% spm_slice_vol.c               2.1
% spm_slice_vol.dll             2.1
% spm_slice_vol.m               2.1
% spm_slice_vol.mexlx           2.1
% spm_slice_vol.mexsol          2.1
% spm_smooth.m                  1.9
% spm_smooth_ui.m               2.4
% spm_sn3d.m                    2.18
% spm_snbasis.m                 2.4
% spm_snbasis_map.m             1.5
% spm_sp.m                      2.8
% spm_spm.m                     2.28
% spm_spm_ui.m                  2.25
% spm_sptop.m                   1.7
% spm_str_manip.m               2.8
% spm_sys_deps.h                2.4
% spm_t2z.m                     2.1
% spm_templates.man             2.5
% spm_transverse.m              2.14
% spm_type.m                    2.3
% spm_u.m                       2.1
% spm_uc.m                      2.3
% spm_unlink.c                  1.2
% spm_unlink.dll                1.2
% spm_unlink.m                  2.1
% spm_unlink.mexlx              1.2
% spm_unlink.mexsol             1.2
% spm_unmap.m                   1.1
% spm_unmap_vol.c               2.1
% spm_unmap_vol.dll             2.1
% spm_unmap_vol.m               2.1
% spm_unmap_vol.mexlx           2.1
% spm_unmap_vol.mexsol          2.1
% spm_vol.m                     2.8
% spm_vol_access.c              2.2
% spm_vol_access.h              2.1
% spm_vol_ecat7.m               2.3
% spm_vol_minc.m                2.7
% spm_vol_utils.c               2.4
% spm_vol_utils.h               2.1
% spm_vrun.m                    2.1
% spm_win32utils.c              2.3
% spm_win32utils.dll            2.3
% spm_win32utils.m              2.1
% spm_write_filtered.m          2.6
% spm_write_plane.m             2.8
% spm_write_sn.m                2.2
% spm_write_vol.m               2.3
% spm_xbrain.m                  2.8
% win32mmap.c                   2.1
% win32mmap.h                   2.1
% 
%                           ----------------
%
% ./apriori
%       brainmask.hdr
%       brainmask.img
%       csf.hdr
%       csf.img
%       gray.hdr
%       gray.img
%       white.hdr
%       white.img
% ./canonical
%       avg152PD.hdr
%       avg152PD.img
%       avg152T1.hdr
%       avg152T1.img
%       avg152T2.hdr
%       avg152T2.img
%       avg305T1.hdr
%       avg305T1.img
%       single_subj_T1.hdr
%       single_subj_T1.img
% ./templates
%       EPI.hdr
%       EPI.img
%       PD.hdr
%       PD.img
%       PET.hdr
%       PET.img
%       T1.hdr
%       T1.img
%       T2.hdr
%       T2.img
%       Transm.hdr
%       Transm.img
%       filT1.hdr
%       filT1.img
% ./rend
%       render_single_subj.mat
%       render_smooth_average.mat
%       render_spm96.mat
% 
%_______________________________________________________________________

help Contents

%=======================================================================
% PROGRAMMERS NOTE:
% This (Contents.m) is the contents file for SPM, used by spm('Ver') to
% recover the version number and copyright information. MatLab's ver
% also uses Contents.m files to identify toolbox versions.
% Line1: Version (first word) & copyright information (rest of line).
% Line2: One line description
%=======================================================================
