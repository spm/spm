%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4899 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.type = 'cfg_entry';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.name = 'Output Variable Name';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.tag = 'name';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.strtype = 's';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.extras = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.num = [1 Inf];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.check = @cfg_check_assignin;
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.help = {'Enter a valid MATLAB variable name. The contents of the input to "Output Item" will be assigned to a variable of this name in MATLAB workspace.'};
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_entry.def = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.type = 'cfg_entry';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.name = 'Output Item';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.tag = 'output';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.strtype = 'e';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.extras = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.num = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.check = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.help = {'The contents that is passed to this input that will be assigned to the workspace variable whose name is given in "Output Variable Name".'};
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_entry.def = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'Pass Output to Workspace';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'cfg_assignin';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'Output Variable Name (cfg_entry)';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).sname = 'Output Item (cfg_entry)';
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @cfg_run_assignin;
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{3}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {
                                                                'Assign a computation result to a workspace variable.'
                                                                'The value entered into "Output Item" will be assigned to a MATLAB workspace variable whose name is specified in "Output Variable Name". If this variable already exists, a new variable name will be generated. This can be useful to assess the results of computations with other MATLAB routines or for debugging.'
                                                                }';
