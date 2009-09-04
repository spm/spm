%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3355 $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_toy{1}.add2{1}.cfg_example_add1.a = '<UNDEFINED>';
matlabbatch{1}.cfg_toy{1}.add2{1}.cfg_example_add1.b = '<UNDEFINED>';
matlabbatch{2}.cfg_toy{1}.cfg_example_div.a = '<UNDEFINED>';
matlabbatch{2}.cfg_toy{1}.cfg_example_div.b(1) = cfg_dep;
matlabbatch{2}.cfg_toy{1}.cfg_example_div.b(1).tname = 'Input b';
matlabbatch{2}.cfg_toy{1}.cfg_example_div.b(1).tgt_spec = {};
matlabbatch{2}.cfg_toy{1}.cfg_example_div.b(1).sname = 'Add1: Add1: a + b';
matlabbatch{2}.cfg_toy{1}.cfg_example_div.b(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.cfg_toy{1}.cfg_example_div.b(1).src_output = substruct('()',{1});
matlabbatch{3}.cfg_basicio.cfg_assignin.name = 'div_sum_ab_c';
matlabbatch{3}.cfg_basicio.cfg_assignin.output(1) = cfg_dep;
matlabbatch{3}.cfg_basicio.cfg_assignin.output(1).tname = 'Output Item';
matlabbatch{3}.cfg_basicio.cfg_assignin.output(1).tgt_spec = {};
matlabbatch{3}.cfg_basicio.cfg_assignin.output(1).sname = 'div: a div b: mod';
matlabbatch{3}.cfg_basicio.cfg_assignin.output(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
matlabbatch{3}.cfg_basicio.cfg_assignin.output(1).src_output = substruct('.','mod');
