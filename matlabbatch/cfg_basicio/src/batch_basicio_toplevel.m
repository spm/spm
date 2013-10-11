%-----------------------------------------------------------------------
% Job saved on 02-Oct-2013 13:46:16 by cfg_util (rev $Rev: 5688 $)
% cfg_basicio BasicIO - Unknown
% menu_cfg ConfGUI - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.type = 'cfg_choice';
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.name = 'BasicIO';
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.tag = 'cfg_basicio';
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.values = {
                                                          '<UNDEFINED>'
                                                          '<UNDEFINED>'
                                                          '<UNDEFINED>'
                                                          }';
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.check = [];
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.rewrite_job = @cfg_basicio_rewrite;
matlabbatch{1}.menu_cfg.menu_struct.conf_choice.help = {'This toolbox contains basic input and output functions. The "Named Input" functions can be used to enter values or file names. These inputs can then be passed on to multiple modules, thereby ensuring all of them use the same input value. Some basic file manipulation is implemented in "Change Directory", "Make Directory", "Move Files". Lists of files can be filtered or splitted into parts using "File Set Filter" and "File Set Split". Output values from other modules can be written out to disk or assigned to MATLAB workspace.'};
