%-----------------------------------------------------------------------
% Job saved on 01-Dec-2015 13:23:09 by cfg_util (rev $Rev: 6627 $)
% spm SPM - SPM12 (12.1)
% cfg_basicio BasicIO - Unknown
% dtijobs DTI tools - Unknown
% impexp_NiftiMrStruct NiftiMrStruct - Unknown
% cfg_logextract LogExt - Unknown
% tbx_dpc dPC Tools - Unknown
% prt PRoNTo - Unknown
% perfusion Perfusion - Unknown
% menu_cfg ConfGUI - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg.menu_entry.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.name = 'File Set';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.tag = 'files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.filter = 'any';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.dir = '';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.num = [0 Inf];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.check = [];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.rewrite_job = [];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.help = {'Select a set of files.'};
matlabbatch{1}.menu_cfg.menu_entry.conf_files.def = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_files.type = 'cfg_files';
matlabbatch{2}.menu_cfg.menu_entry.conf_files.name = 'Output directory';
matlabbatch{2}.menu_cfg.menu_entry.conf_files.tag = 'outdir';
matlabbatch{2}.menu_cfg.menu_entry.conf_files.filter = 'dir';
matlabbatch{2}.menu_cfg.menu_entry.conf_files.ufilter = '.*';
matlabbatch{2}.menu_cfg.menu_entry.conf_files.dir = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_files.num = [0 1];
matlabbatch{2}.menu_cfg.menu_entry.conf_files.check = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_files.rewrite_job = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_files.help = {'Output files will be placed in this folder. Leave empty to put them into the same folder as the original files.'};
matlabbatch{2}.menu_cfg.menu_entry.conf_files.def = [];
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.type = 'cfg_menu';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.name = 'Keep original files';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.tag = 'keep';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.labels = {
                                                       'Yes'
                                                       'No'
                                                       }';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.values = {
                                                       true
                                                       false
                                                       }';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.check = [];
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.rewrite_job = [];
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.help = {};
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.def = [];
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.name = 'Gzip Files';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.tag = 'cfg_gzip_files';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1) = cfg_dep('Files: File Set (cfg_files)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}));
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1) = cfg_dep('Files: Output directory (cfg_files)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}));
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1) = cfg_dep('Menu: Keep original files (cfg_menu)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}));
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.prog = @cfg_run_gzip_files;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.vout = @cfg_vout_gzip_files;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.check = [];
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.rewrite_job = [];
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.help = {'Gzip each file in a set of files.'};
