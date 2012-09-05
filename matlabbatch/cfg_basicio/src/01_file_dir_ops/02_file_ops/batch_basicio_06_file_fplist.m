%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4899 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg.menu_entry.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.name = 'Directory';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.tag = 'dir';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.filter = 'dir';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.dir = '';
matlabbatch{1}.menu_cfg.menu_entry.conf_files.num = [1 Inf];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.check = [];
matlabbatch{1}.menu_cfg.menu_entry.conf_files.help = {'Directory to select files from.'};
matlabbatch{1}.menu_cfg.menu_entry.conf_files.def = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.type = 'cfg_entry';
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.name = 'Filter';
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.tag = 'filter';
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.strtype = 's';
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.extras = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.num = [1 Inf];
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.check = [];
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.help = {'A regular expression to filter files.'};
matlabbatch{2}.menu_cfg.menu_entry.conf_entry.def = [];
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.type = 'cfg_menu';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.name = 'Descend into subdirectories';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.tag = 'rec';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.labels = {
                                                       'Yes'
                                                       'No'
                                                       }';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.values = {
                                                       'FPListRec'
                                                       'FPList'
                                                       }';
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.check = [];
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.help = {'Files can be selected from the specified directory only or from the specified directory and all its subdirectories.'};
matlabbatch{3}.menu_cfg.menu_entry.conf_menu.def = [];
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.name = 'File Selector (Batch Mode)';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.tag = 'file_fplist';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1).sname = 'Directory (cfg_files)';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1).sname = 'Filter (cfg_entry)';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1) = cfg_dep;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1).sname = 'Menu: Descend into subdirectories (cfg_menu)';
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.val{3}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.prog = @cfg_run_file_fplist;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.vout = @cfg_vout_file_fplist;
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.check = [];
matlabbatch{4}.menu_cfg.menu_struct.conf_exbranch.help = {'Select files from a directory using cfg_getfile(''FPList'',...).'};
