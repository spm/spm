%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3921 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.type = 'cfg_choice';
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.name = 'BasicIO';
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.tag = 'cfg_basicio';
%%
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.values = {
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
                                                          '<UNDEFINED>'
                                                          '<UNDEFINED>'
                                                          '<UNDEFINED>'
                                                          }';
%%
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.check = [];
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.help = {'This toolbox contains basic input and output functions. The "Named Input" functions can be used to enter values or file names. These inputs can then be passed on to multiple modules, thereby ensuring all of them use the same input value. Some basic file manipulation is implemented in "Change Directory", "Make Directory", "Move Files". Lists of files can be filtered or splitted into parts using "File Set Filter" and "File Set Split". Output values from other modules can be written out to disk or assigned to MATLAB workspace.'};
matlabbatch{2}.menu_cfg.gencode_gen.gencode_fname = 'cfg_cfg_basicio.m';
matlabbatch{2}.menu_cfg.gencode_gen.gencode_dir = '<UNDEFINED>';
matlabbatch{2}.menu_cfg.gencode_gen.gencode_var(1) = cfg_dep;
matlabbatch{2}.menu_cfg.gencode_gen.gencode_var(1).tname = 'Root node of config';
matlabbatch{2}.menu_cfg.gencode_gen.gencode_var(1).tgt_spec = {};
matlabbatch{2}.menu_cfg.gencode_gen.gencode_var(1).sname = 'BasicIO (cfg_repeat)';
matlabbatch{2}.menu_cfg.gencode_gen.gencode_var(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.menu_cfg.gencode_gen.gencode_var(1).src_output = substruct('()',{1});
matlabbatch{2}.menu_cfg.gencode_gen.gencode_opts.gencode_o_def = true;
matlabbatch{2}.menu_cfg.gencode_gen.gencode_opts.gencode_o_mlb = false;
matlabbatch{2}.menu_cfg.gencode_gen.gencode_opts.gencode_o_path = false;
