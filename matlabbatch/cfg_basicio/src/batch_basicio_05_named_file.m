%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 2787 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.type = 'cfg_entry';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.name = 'Input Name';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.tag = 'name';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.strtype = 's';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.extras = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.num = [1 Inf];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.check = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.help = {'Enter a name for this file selection. This name will be displayed in the ''Dependency'' listing as output name.'};
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.def = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.type = 'cfg_files';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.name = 'File Set';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.tag = 'files';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.filter = 'any';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.ufilter = '.*';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.dir = '';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.num = [0 Inf];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.check = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.help = {'Select a set of files.'};
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_files.def = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.type = 'cfg_repeat';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.name = 'File Sets';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.tag = 'files';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1) = cfg_dep;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tname = 'Values Item';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tgt_spec = {};
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).sname = 'File selection (cfg_files)';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_output = substruct('()',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.num = [1 Inf];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.forcestruct = false;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.check = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.help = {'Select one or more sets of files. Each set can be passed separately as a dependency to other modules.'};
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'Named File Selector';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'cfg_named_file';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'Input Name (cfg_entry)';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).sname = 'File Sets (cfg_repeat)';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @cfg_run_named_file;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = @cfg_vout_named_file;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {
                                                                'Named File Selector allows to select sets of files that can be referenced as common input by other modules.'
                                                                'In addition to file outputs, an index vector is provided that can be used to separate the files again using ''File Set Split'' module. This can be useful if the same sets of files have to be processed separately in one module and as a single set in another.'
                                                                }';
