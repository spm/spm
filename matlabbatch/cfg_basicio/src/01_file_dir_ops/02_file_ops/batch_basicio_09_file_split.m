%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4899 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.type = 'cfg_entry';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.name = 'File Set Name';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.tag = 'name';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.strtype = 's';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.extras = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.num = [1 Inf];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.check = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.help = {'Enter a name for this file selection. This name will be displayed in the ''Dependency'' listing as output name.'};
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.def = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.type = 'cfg_entry';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.name = 'Selection Index';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.tag = 'index';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.strtype = 'n';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.extras = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.num = [1 Inf];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.check = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.help = {
                                                            'Enter the index vector of files belonging to this set. E.g. to select files 1 to 10 from the input file set, enter [1:10].'
                                                            'To split the combined list of all files selected by a "Named File Selector", enter the index vectors here as dependency.'
                                                            }';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.def = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.type = 'cfg_repeat';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.name = 'Output File Sets';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.tag = 'filesets';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1) = cfg_dep;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tname = 'Values Item';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tgt_spec = {};
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).sname = 'Selection Index (cfg_entry)';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_output = substruct('()',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.num = [1 Inf];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.forcestruct = false;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.check = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_repeat.help = {'For each file set to be created, enter an index vector that selects files from the input list. An additional output file set will be created matching files not selected by any index.'};
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.type = 'cfg_files';
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.name = 'Input File Set';
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.tag = 'files';
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.filter = 'any';
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.ufilter = '.*';
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.dir = '';
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.num = [1 Inf];
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.check = [];
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.help = {'The input file set will be split at the indices given in the ''#files per set'' collection.'};
matlabbatch{4}.menu_cfg{1}.menu_entry{1}.conf_files.def = [];
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'File Set Split';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'cfg_file_split';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'File Set Name (cfg_entry)';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).sname = 'Input File Set (cfg_files)';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{3}(1) = cfg_dep;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{3}(1).tname = 'Val Item';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{3}(1).tgt_spec = {};
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{3}(1).sname = 'File Sets (cfg_repeat)';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{3}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{3}(1).src_output = substruct('()',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @cfg_run_file_split;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = @cfg_vout_file_split;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {
                                                                'Split a list of files into multiple parts.'
                                                                'With this utility, a list of files can be split into parts based on index vectors. Each index vector is a list of numbers where number N selects the Nth file from the input list. Index vectors may overlap, so files may belong to more than one part. Files not matching any index will be returned in a separate list. If an index is larger than the number of files passed it will be ignored.'
                                                                }';
