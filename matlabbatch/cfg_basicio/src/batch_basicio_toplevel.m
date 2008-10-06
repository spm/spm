%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 2305 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.type = 'cfg_choice';
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.name = 'BasicIO';
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.tag = 'cfg_basicio';
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.values = {
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
                                                                '<UNDEFINED>'
}';
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.num = double([0 Inf]);
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.forcestruct = logical(true);
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.check = double([]);
matlabbatch{1}.menu_cfg{1}.menu_struct{1}.conf_choice.help = {'This toolbox contains basic input and output functions. The "Named Input" functions can be used to enter values or file names. These inputs can then be passed on to multiple modules, thereby ensuring all of them use the same input value. Some basic file manipulation is implemented in "Change Directory", "Make Directory", "Move Files". Lists of files can be filtered or splitted into parts using "File Set Filter" and "File Set Split". Output values from other modules can be written out to disk or assigned to MATLAB workspace.'};
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_fname = 'cfg_cfg_basicio.m';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_dir = {'/export/spm-devel/matlabbatch/trunk/cfg_basicio/'};
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1) = cfg_dep;
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).tname = 'Root node of config';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).tgt_spec = {};
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).sname = 'BasicIO (cfg_repeat)';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(1).type = '.';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(1).subs = 'val';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(2).type = '{}';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(2).subs{1} = double(1);
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(3).type = '.';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(3).subs = 'val';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(4).type = '{}';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(5).type = '.';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(5).subs = 'val';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(6).type = '{}';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_output(1).type = '()';
matlabbatch{2}.menu_cfg{1}.gencode_gen.gencode_var(1).src_output(1).subs{1} = double(1);
