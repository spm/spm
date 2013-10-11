%-----------------------------------------------------------------------
% Job saved on 25-Sep-2013 15:22:57 by cfg_util (rev $Rev: 5685 $)
% spm SPM - SPM12b (beta)
% cfg_basicio BasicIO - Unknown
% cfg_logextract LogExt - Unknown
% tbx_dpc dPC Tools - Unknown
% perfusion Perfusion - Unknown
% dtijobs DTI tools - Unknown
% impexp_NiftiMrStruct NiftiMrStruct - Unknown
% prt PRoNTo - Unknown
% menu_cfg ConfGUI - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg.menu_entry.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.name = 'File Set';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.tag = 'files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.filter = '\.gz$';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.dir = '';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.num = [0 Inf];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.check = [];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.help = {'Select a set of files.'};
matlabbatch{1}.menu_cfg.menu_entry.conf_files.def = [];
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.name = 'GunZip Files';
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.tag = 'cfg_gunzip_files';
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.val{1}(1) = cfg_dep('Files: File Set (cfg_files)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}));
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.prog = @(job)gunzip(job.files);
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.vout = @cfg_vout_gunzip_files;
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.check = [];
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.help = {'GunZip each file in a set of files.'};
