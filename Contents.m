% SPM98a (c) 1991,1994-1996,1998: The Wellcome Department of Cognitive Neurology
% Statistical Parametric Mapping - SPM'98 developers alpha release (ML5.2)
%_______________________________________________________________________
% 
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  Institute of Neurology, University College London
%
% John Ashburner, Karl Friston, Andrew Holmes, Jean-Baptiste Poline
%_______________________________________________________________________
%
% This Contents.m file holds the version ID for this release of Matlab,
% and contains a manifest of the included functions and their version numbers.
%
% SPM98 is written for Matlab v5.2 under UNIX
% (compiled binaries of the external functions are provided for Solaris2 only)
%
% See spm.man for details of this release.
% See the README for information on installation and getting started.
% See spm_motd.man for last minute release details.
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%
%
% SPM98a - Manifest
%-----------------------------------------------------------------------
% Contents.m			2.1	
% Grid.mat			1.1
% MIP.mat			1.2
% README			2.1
% Split.mat			1.1
% connex.h			1.2
% dbh.h				(Mayo clinic database sub-definitions)
% render.mat			1.1	
% spm.m				2.19
% spm.man			2.1
% spm_AnCova.m			2.3
% spm_Bcdf.m			2.1
% spm_Bpdf.m			2.1
% spm_DesMtx.m			2.5
% spm_DesMtxSca.m		2.1
% spm_DesRep.m			2.4
% spm_F.m			1.2
% spm_F.man			1.2
% spm_Fcdf.m			2.1
% spm_FcnTpl.txt		1.1
% spm_Fpdf.m			2.1
% spm_Gcdf.m			2.1
% spm_Gpdf.m			2.1
% spm_Icdf.m			2.1
% spm_Ipdf.m			2.1
% spm_MAKE			2.4
% spm_Ncdf.m			2.1
% spm_Npdf.m			2.1
% spm_P.m			1.2
% spm_Pcdf.m			2.1
% spm_Pec.m			1.1
% spm_Pec_resels.m		1.1
% spm_Pkn.m			1.3
% spm_Pn.m			1.1
% spm_Ppdf.m			2.1
% spm_Pz.m			1.1
% spm_RandFX.man		1.3
% spm_SpUtil.m			2.8
% spm_Tcdf.m			2.1
% spm_Tpdf.m			2.1
% spm_VOI.m			1.1
% spm_Volt.m			1.1
% spm_Volt_W.m			1.2
% spm_W.m			1.2
% spm_XYZreg.m			2.1
% spm_XYZreg_Ex1.m		2.1
% spm_XYZreg_Ex2.m		2.1
% spm_Xcdf.m			2.1
% spm_Xpdf.m			2.1
% spm_add.c			2.5
% spm_add.m			2.6
% spm_add.mexsol		2.5
% spm_adjmean_fmri_ui.m		2.3
% spm_adjmean_ui.m		2.3
% spm_affsub3.m			2.2
% spm_append.m			2.1
% spm_atranspa.c		1.2
% spm_atranspa.m		1.2
% spm_atranspa.mexsol		1.2
% spm_brainwarp.c		2.1
% spm_brainwarp.m		2.1
% spm_brainwarp.mexsol		2.1
% spm_check_registration.m	1.1
% spm_chi2_plot.m		1.2
% spm_choose.m			2.3
% spm_clf.m			2.1
% spm_clusters.c		1.1
% spm_clusters.m		1.1
% spm_clusters.mexsol		1.1
% spm_conv.m			1.4
% spm_conv_vol.c		1.7
% spm_conv_vol.m		1.3
% spm_conv_vol.mexsol		1.7
% spm_coregister.m		2.1
% spm_create_image.m		2.3
% spm_dbm.m			2.1
% spm_dctmtx.m			1.3
% spm_defaults.m		2.3
% spm_defaults_edit.m		2.1
% spm_detrend.m			1.2
% spm_display.m			2.2
% spm_dummy.m			2.1
% spm_efmri.man			1.2
% spm_en.m			1.1
% spm_extract.m			2.2
% spm_figure.m			2.6
% spm_fmri.man			1.2
% spm_fmri_spm_ui.m		1.35
% spm_format.man		1.2
% spm_fzero.m			2.1
% spm_fzeroIB.m			2.1
% spm_get.m			2.14
% spm_get_space.m		2.1
% spm_global.c			2.2
% spm_global.m			2.1
% spm_global.mexsol		2.2
% spm_graph.m			1.18
% spm_graph_ui.m		2.1
% spm_grid.m			1.1
% spm_header_edit.m		1.2
% spm_help.m			2.12
% spm_hread.m			2.1
% spm_hrf.m			1.4
% spm_hwrite.m			1.4
% spm_image.m			1.4
% spm_image.man			1.1
% spm_imatrix.m			1.2
% spm_imcalc.m			2.6
% spm_imcalc_ui.m		2.6
% spm_input.m			2.23
% spm_invBcdf.m			2.1
% spm_invFcdf.m			2.1
% spm_invGcdf.m			2.1
% spm_invIcdf.m			2.1
% spm_invNcdf.m			2.1
% spm_invPcdf.m			2.1
% spm_invTcdf.m			2.1
% spm_invXcdf.m			2.1
% spm_invkcdf.m			1.1
% spm_k.m			1.1
% spm_kcdf.m			1.1
% spm_kpdf.m			1.1
% spm_kronutil.c		1.1
% spm_kronutil.m		2.1
% spm_kronutil.mexsol		1.1
% spm_lambda.m			1.7
% spm_list_files.c		2.1
% spm_list_files.m		1.2
% spm_list_files.mexsol		2.1
% spm_load.m			1.3
% spm_log.m			2.2
% spm_make_lookup.c		2.1
% spm_map.h			2.1
% spm_map.m			1.5
% spm_map_vol.c			1.9
% spm_map_vol.m			1.1
% spm_map_vol.mexsol		1.9
% spm_mapping.c			2.1
% spm_mask.m			2.5
% spm_mat.man			1.1
% spm_matrix.m			1.1
% spm_matx.m			2.4
% spm_max.c			1.1
% spm_max.m			1.1
% spm_max.mexsol		1.1
% spm_maxima.m			1.9
% spm_mean_ui.m			2.4
% spm_meanby.m			2.1
% spm_methods.man		1.1
% spm_min_Pn.m			1.1
% spm_min_Pz.m			1.1
% spm_mip.m			1.10
% spm_mip_ui.m			2.1
% spm_modality.man		1.1
% spm_motd.man			2.1
% spm_mvNpdf.m			2.1
% spm_nCr.m			2.1
% spm_orthviews.m		2.2
% spm_pF.m			2.1
% spm_pet.man			1.2
% spm_picture.m			2.1
% spm_print.m			1.3
% spm_progress_bar.m		1.4
% spm_project.c			1.6
% spm_project.m			1.1
% spm_project.mexsol		1.6
% spm_projections.m		1.15
% spm_projections.man		1.1
% spm_projectionsF_ui.m		2.1
% spm_projections_ui.m		2.1
% spm_readXA.m			1.2
% spm_read_vols.m		2.4
% spm_realign.m			2.4
% spm_realign.man		1.3
% spm_render.m			2.1
% spm_render_vol.c		1.5
% spm_render_vol.m		1.1
% spm_render_vol.mexsol		1.5
% spm_renviews.m		1.4
% spm_resize.m			1.1
% spm_resss.m			2.4
% spm_results.m			2.1
% spm_results_ui.m		2.2
% spm_sample_vol.c		1.11
% spm_sample_vol.m		1.2
% spm_sample_vol.mexsol		1.11
% spm_sections.m		2.2
% spm_segment.m			2.2
% spm_slice_vol.c		1.8
% spm_slice_vol.m		1.2
% spm_slice_vol.mexsol		1.8
% spm_smooth.m			1.9
% spm_smooth.man		1.1
% spm_smooth_ui.m		1.3
% spm_sn3d.m			2.4
% spm_sn3d.man			1.2
% spm_snbasis.m			2.1
% spm_snbasis_map.m		1.5
% spm_sp.m			2.3
% spm_spm.m			2.4
% spm_spm.man			1.1
% spm_spm_ui.m			2.14
% spm_sptop.m			1.7
% spm_str_manip.m		2.3
% spm_svd.m			1.4
% spm_svd.man			1.1
% spm_svd_ui.m			1.3
% spm_t2z.m			2.1
% spm_transverse.m		1.10
% spm_type.m			2.2
% spm_unlink.c			1.2
% spm_unlink.m			1.2
% spm_unlink.mexsol		1.2
% spm_unmap.m			1.1
% spm_unmap_vol.c		1.4
% spm_unmap_vol.m		1.1
% spm_unmap_vol.mexsol		1.4
% spm_vol.m			2.4
% spm_vol_ecat7.m		1.1
% spm_vol_minc.m		2.4
% spm_vol_utils.c		2.1
% spm_vol_utils.h		2.1
% spm_vol_utils2.c		2.1
% spm_write.m			2.2
% spm_write_filtered.m		1.1
% spm_write_plane.m		2.2
% spm_write_sn.m		1.10
% spm_write_vol.m		2.2
% spm_z.m			1.1
% templates.man			2.1
% volume.h			1.7
% 
%                           ----------------
%
% ./apriori
% 	brainmask.hdr
% 	brainmask.img
% 	csf.hdr
% 	csf.img
% 	gray.hdr
% 	gray.img
% 	white.hdr
% 	white.img
% ./templates
% 	EPI.hdr
% 	EPI.img
% 	PET.hdr
% 	PET.img
% 	T1.hdr
% 	T1.img
% 	T2.hdr
% 	T2.img
% 	Transm.hdr
% 	Transm.img
% ./coreg
% 	EPI.hdr
% 	EPI.img
% 	PET.hdr
% 	PET.img
% 	T1.hdr
% 	T1.img
% 	T2.hdr
% 	T2.img
% 	Transm.hdr
% 	Transm.img
%_______________________________________________________________________



%=======================================================================
% PROGRAMMERS NOTE:
% This (Contents.m) is the contents file for SPM, used by spm('Ver') to
% recover the version number and copyright information. MatLab's ver
% also uses Contents.m files to identify toolbox versions.
% Line1: Version (first word) & copyright information (rest of line).
% Line2: One line description
%=======================================================================
