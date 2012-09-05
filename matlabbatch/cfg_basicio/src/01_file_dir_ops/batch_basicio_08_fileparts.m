%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4899 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg.menu_entry.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.name = 'Files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.tag = 'files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.filter = 'any';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.dir = [];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.num = [1 Inf];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.check = [];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.help = {'Enter file names.'};
matlabbatch{1}.menu_cfg.menu_entry.conf_files.def = [];
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.name = 'Get Pathnames';
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.tag = 'cfg_fileparts';
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.val{1}(1) = cfg_dep('Files: Files (cfg_files)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}));
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.prog = @cfg_run_fileparts;
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.vout = @cfg_vout_fileparts;
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.check = [];
matlabbatch{2}.menu_cfg.menu_struct.conf_exbranch.help = {'Split a list of path names into lists of parent directory, name and extension. The list of parent directories can be passed to any other file selector. The list of names and extensions can only be passed to evaluated inputs, not to string inputs.'};
