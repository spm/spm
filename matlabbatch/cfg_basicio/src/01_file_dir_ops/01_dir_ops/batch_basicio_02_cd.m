%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4899 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.name = 'Directory';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.tag = 'dir';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.filter = 'dir';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.dir = '';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.num = [1 1];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.check = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.help = {'New working directory.'};
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.def = [];
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'Change Directory';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'cfg_cd';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'Directory (cfg_files)';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @cfg_run_cd;
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = [];
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {'Change working directory.'};
