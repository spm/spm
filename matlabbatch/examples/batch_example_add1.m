%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 1184 $)
%-----------------------------------------------------------------------
matlabbatch{1, 1}.cfg_toy{1, 1}.add2{1, 1}.cfg_example_add1.a = double(1);
matlabbatch{1, 1}.cfg_toy{1, 1}.add2{1, 1}.cfg_example_add1.b = double(2);
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.name = 'add1';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1) = cfg_dep;
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).tname = 'Output Item';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).sname = 'Add1: a + b';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 1).type = '.';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 1).subs = 'val';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 2).type = '{}';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 2).subs{1, 1} = double(1);
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 3).type = '.';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 3).subs = 'val';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 4).type = '{}';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 4).subs{1, 1} = double(1);
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 5).type = '.';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 5).subs = 'val';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 6).type = '{}';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_exbranch(1, 6).subs{1, 1} = double(1);
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_output(1).type = '()';
matlabbatch{1, 2}.cfg_basicio{1, 1}.cfg_assignin.output(1).src_output(1).subs{1, 1} = double(1);
