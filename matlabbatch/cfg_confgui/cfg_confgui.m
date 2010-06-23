function menu_cfg = cfg_confgui

% This function describes the user defined fields for each kind of
% cfg_item and their layout in terms of cfg_items. Thus, the
% configuration system can be used to generate code for new configuration
% files itself.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_confgui.m 3944 2010-06-23 08:53:40Z volkmar $

rev = '$Rev: 3944 $'; %#ok

%% Declaration of fields

% Name
%-----------------------------------------------------------------------
conf_name         = cfg_entry;
conf_name.name    = 'Name';
conf_name.tag     = 'name';
conf_name.strtype = 's';
conf_name.num     = [1 Inf];
conf_name.help    = {'Display name of configuration item.'};

% Tag
%-----------------------------------------------------------------------
conf_tag         = cfg_entry;
conf_tag.name    = 'Tag';
conf_tag.tag     = 'tag';
conf_tag.strtype = 's';
conf_tag.num     = [1 Inf];
conf_tag.help    = {['Tag of configuration item.', 'This will be used as tag' ...
                    ' of the generated output and (depending on item class)' ...
                    ' appear in the input structure to the computation' ...
                    ' function.']};

% Val Item
%-----------------------------------------------------------------------
conf_val_item         = cfg_entry;
conf_val_item.name    = 'Val Item';
conf_val_item.tag     = 'val';
conf_val_item.strtype = 'e';
conf_val_item.num     = [1 1];
conf_val_item.help    = {'Val of configuration item.', 'This should be a dependency to another cfg_item object.'};

% Val
%-----------------------------------------------------------------------
conf_val         = cfg_repeat;
conf_val.name    = 'Val';
conf_val.tag     = 'val';
conf_val.values  = {conf_val_item};
conf_val.num     = [1 Inf];
conf_val.help    = {'Val of configuration item.', 'A collection of cfg_item objects to be assembled in a cfg_(ex)branch. Each item in this list needs to have a unique tag.'};

% Val with exactly one element 
%-----------------------------------------------------------------------
conf_val_single = conf_val;
conf_val_single.val = {conf_val_item};
conf_val_single.num = [1 1];

% Check
%-----------------------------------------------------------------------
conf_check         = cfg_entry;
conf_check.name    = 'Check';
conf_check.tag     = 'check';
conf_check.strtype = 'f';
conf_check.num     = [0 Inf];
conf_check.help    = {'Check function (handle).', ...
    ['This function will be called during all_set, before a job can be run. ', ...
    'It receives the harvested configuration tree rooted at the current item as input. ', ...
    'If the input is ok, it should return an empty string. Otherwise, ', ...
    'its output should be a string that describes why input is not correct or consistent.'], ...
    ['Note that the check function will be called only if all dependencies are resolved. ', ...
    'This will usually be at the time just before the job is actually run.']};

% Help paragraph
%-----------------------------------------------------------------------
conf_help_par         = cfg_entry;
conf_help_par.name    = 'Paragraph';
conf_help_par.val     = {''};
conf_help_par.tag     = 'help';
conf_help_par.strtype = 's';
conf_help_par.num     = [0 Inf];
conf_help_par.help    = {'Help paragraph.', 'Enter a string which will be formatted as a separate paragraph.'};

% Help
%-----------------------------------------------------------------------
conf_help         = cfg_repeat;
conf_help.name    = 'Help';
conf_help.tag     = 'help';
conf_help.values  = {conf_help_par};
conf_help.num     = [0 Inf];
conf_help.help    = {'Help text.', 'Each help text consists of a number of paragraphs.'};

% Def
%-----------------------------------------------------------------------
conf_def         = cfg_entry;
conf_def.name    = 'Def';
conf_def.tag     = 'def';
conf_def.strtype = 'f';
conf_def.num     = [0 1];
conf_def.help    = {'Default settings for configuration item.', ...
         ['This should be a function handle to a function accepting '...
          'both zero and one free argument. It will be called by setval ' ...
          'with zero free arguments to retrieve a default value and with ' ...
          'one argument to set a default.'], ...
                    ['If the default function has a registry like call ' ...
                    'syntax '], ...
                    'defval = get_defaults(''some.key'')', ...
                    'get_defaults(''some.key'', defval)', ...
                    'then the function handle should look like this', ...
                    '@(defval)get_defaults(''some.key'', defval{:})', ...
                    ['Matlabbatch will wrap the second argument in a cell ' ...
                    'when calling this handle, thereby effectively passing ' ...
                    'no second argument to retrieve a value and passing ' ...
                    'the default value when setting it.']};

% Hidden
%-----------------------------------------------------------------------
conf_hidden         = cfg_menu;
conf_hidden.name    = 'Hidden';
conf_hidden.tag     = 'hidden';
conf_hidden.labels  = {'False', 'True'};
conf_hidden.values  = {false, true};
conf_hidden.help    = {'''Hidden'' status of configuration item.', 'If you want to hide an item to the user, set this to true. However, the item will be harvested, and a job can only run if this item is all_set. To hide an item is mostly useful for constants, that do not need to be changed by the user.'};

% Forcestruct
%-----------------------------------------------------------------------
conf_forcestruct         = cfg_menu;
conf_forcestruct.name    = 'Forcestruct';
conf_forcestruct.tag     = 'forcestruct';
conf_forcestruct.labels  = {'False', 'True'};
conf_forcestruct.values  = {false, true};
conf_forcestruct.help    = {'Forcestruct flag.', 'Sometimes the default harvest behaviour of a cfg_repeat object is not what one wants to have: if there is only one repeatable item, it returns a cell and discards its own tag. If one adds a new configuration to this repeat, then suddenly the harvested cfg_repeat has the tag of the repeat item and a cell with struct arrays as members. If you want to force this latter behaviour also if there is only one item in the ''values'' field, then set this flag to ''true''.'};

% Values Item
%-----------------------------------------------------------------------
conf_values_item         = cfg_entry;
conf_values_item.name    = 'Values Item';
conf_values_item.tag     = 'values';
conf_values_item.strtype = 'e';
conf_values_item.num     = [1 Inf];
conf_values_item.help    = {'Value of configuration item.', 'For cfg_menus, this is an arbitrary value, for cfg_repeat/cfg_choice items, this should be a dependency to another cfg_item object.'};

% Values
%-----------------------------------------------------------------------
conf_values         = cfg_repeat;
conf_values.name    = 'Values';
conf_values.tag     = 'values';
conf_values.values  = {conf_values_item};
conf_values.num     = [1 Inf];
conf_values.help    = {'Values of configuration item. If this is a cfg_repeat/cfg_choice and there is more than one item in this list, each item needs to have a unique tag.'};

% Label Item
%-----------------------------------------------------------------------
conf_labels_item         = cfg_entry;
conf_labels_item.name    = 'Label';
conf_labels_item.tag     = 'labels';
conf_labels_item.strtype = 's';
conf_labels_item.num     = [1 Inf];
conf_labels_item.help    = {'Label of menu item.', 'This is a string which will become a label entry.'};

% Labels
%-----------------------------------------------------------------------
conf_labels         = cfg_repeat;
conf_labels.name    = 'Labels';
conf_labels.tag     = 'labels';
conf_labels.values  = {conf_labels_item};
conf_labels.num     = [1 Inf];
conf_labels.help    = {'Labels of configuration item.', 'This is a collection of strings - each string will become a label entry.'};

% Filter
%-----------------------------------------------------------------------
conf_filter         = cfg_entry;
conf_filter.name    = 'Filter';
conf_filter.tag     = 'filter';
conf_filter.strtype = 's'; 
conf_filter.num     = [0 Inf];
conf_filter.help    = {'Filter for files.', ...
    ['This filter will not display in the file selector filter field, '...
    'it will be used prior to displaying files. The following special '...
    'types are supported:'], ...
    '* ''any''', ...
    '* ''batch''', ...
    '* ''dir''', ...
    '* ''image''', ...
    '* ''mat''', ...
    '* ''mesh''', ...
    '* ''nifti''', ...
    '* ''xml''', ...
    '* any regular expression filter.'};

% Ufilter
%-----------------------------------------------------------------------
conf_ufilter         = cfg_entry;
conf_ufilter.name    = 'Ufilter';
conf_ufilter.tag     = 'ufilter';
conf_ufilter.strtype = 's';
conf_ufilter.num     = [0 Inf];
conf_ufilter.help    = {'Filter for files.', 'This filter is a regexp to filter filenames that survive .filter field.'};

% Dir
%-----------------------------------------------------------------------
conf_dir         = cfg_entry;
conf_dir.name    = 'Dir';
conf_dir.tag     = 'dir';
conf_dir.strtype = 'e'; % TODO: This should really be a cell string type
conf_dir.num     = [0 Inf];
conf_dir.help    = {'Dir field.', 'Initial directory for file selector.'};

% Num
%-----------------------------------------------------------------------
conf_num_any         = cfg_entry;
conf_num_any.name    = 'Num';
conf_num_any.tag     = 'num';
conf_num_any.strtype = 'w';
conf_num_any.num     = [0 Inf];
conf_num_any.help    = {'Num field.', ...
                    ['Specify how many dimensions this ' ...
                    'item must have and the size of each dimension. An ' ...
                    'empty num field means no restriction on number of ' ...
                    'dimensions. Inf in any dimension means no upper limit ' ...
                    'on size of this dimension.'], ...
                    ['Note that for strtype ''s'' inputs, num is interpreted ' ...
                    'as a 2-vector [min max], allowing for a 1-by-min...max ' ...
                    'string to be entered.']};

% Num
%-----------------------------------------------------------------------
conf_num         = cfg_entry;
conf_num.name    = 'Num';
conf_num.tag     = 'num';
conf_num.strtype = 'w';
conf_num.num     = [1 2];
conf_num.help    = {'Num field.', 'Specify how many items (min and max) must be entered to yield a valid input. min = 0 means zero or more, max = Inf means no upper limit on number.'};

% Strtype
%-----------------------------------------------------------------------
conf_strtype        = cfg_menu;
conf_strtype.name   = 'Strtype';
conf_strtype.tag    = 'strtype';
conf_strtype.labels = {'String (s)', ...
                    'Evaluated (e)', ...
                    'Natural number (1..n) (n)', ...
                    'Whole number (0..n) (w)', ...
                    'Integer (i)', ...
                    'Real number (r)', ...
                    'Function handle (f)', ...
                    'Condition vector (c)', ...
                    'Contrast matrix (x)', ...
                    'Permutation (p)'}; 
conf_strtype.values = {'s','e','n','w','i','r','f','c','x','p'};
conf_strtype.help   = {'Strtype field.', 'This type describes how an evaluated input should be treated. Type checking against this type will be performed during subscript assignment.'};

% Extras
%-----------------------------------------------------------------------
conf_extras         = cfg_entry;
conf_extras.name    = 'Extras';
conf_extras.tag     = 'extras';
conf_extras.strtype = 'e';
conf_extras.num     = [0 Inf];
conf_extras.help    = {'Extras field.', 'Extra information that may be used to evaluate ''strtype''.'};

% Prog
%-----------------------------------------------------------------------
conf_prog         = cfg_entry;
conf_prog.name    = 'Prog';
conf_prog.tag     = 'prog';
conf_prog.strtype = 'f';
conf_prog.num     = [1 Inf];
conf_prog.help    = {'Prog function (handle).', 'This function will be called to run a job. It receives the harvested configuration tree rooted at the current item as input. If it produces output, this should be a single variable. This variable can be a struct, cell or whatever is appropriate. To pass references to it a ''vout'' function has to be implemented that describes the virtual outputs.'};

% Vout
%-----------------------------------------------------------------------
conf_vout         = cfg_entry;
conf_vout.name    = 'Vout';
conf_vout.tag     = 'vout';
conf_vout.strtype = 'f';
conf_vout.num     = [0 Inf];
conf_vout.help    = {'Vout function (handle).', 'This function will be called during harvest, if all inputs to a job are set. It receives the harvested configuration tree rooted at the current item as input. Its output should be an array of cfg_dep objects, containing subscript indices into the output variable that would result when running this job. Note that no dependencies are resolved here.'};

%% Declaration of item classes

% Branch
%-----------------------------------------------------------------------
conf_class_branch        = cfg_const;
conf_class_branch.name   = 'Branch';
conf_class_branch.tag    = 'type';
conf_class_branch.val    = {'cfg_branch'};
conf_class_branch.hidden = true;
conf_class_branch.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Choice
%-----------------------------------------------------------------------
conf_class_choice        = cfg_const;
conf_class_choice.name   = 'Choice';
conf_class_choice.tag    = 'type';
conf_class_choice.val    = {'cfg_choice'};
conf_class_choice.hidden = true;
conf_class_choice.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Const
%-----------------------------------------------------------------------
conf_class_const        = cfg_const;
conf_class_const.name   = 'Const';
conf_class_const.tag    = 'type';
conf_class_const.val    = {'cfg_const'};
conf_class_const.hidden = true;
conf_class_const.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Entry
%-----------------------------------------------------------------------
conf_class_entry        = cfg_const;
conf_class_entry.name   = 'Entry';
conf_class_entry.tag    = 'type';
conf_class_entry.val    = {'cfg_entry'};
conf_class_entry.hidden = true;
conf_class_entry.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Exbranch
%-----------------------------------------------------------------------
conf_class_exbranch        = cfg_const;
conf_class_exbranch.name   = 'Exbranch';
conf_class_exbranch.tag    = 'type';
conf_class_exbranch.val    = {'cfg_exbranch'};
conf_class_exbranch.hidden = true;
conf_class_exbranch.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Files
%-----------------------------------------------------------------------
conf_class_files        = cfg_const;
conf_class_files.name   = 'Files';
conf_class_files.tag    = 'type';
conf_class_files.val    = {'cfg_files'};
conf_class_files.hidden = true;
conf_class_files.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Menu
%-----------------------------------------------------------------------
conf_class_menu        = cfg_const;
conf_class_menu.name   = 'Menu';
conf_class_menu.tag    = 'type';
conf_class_menu.val    = {'cfg_menu'};
conf_class_menu.hidden = true;
conf_class_menu.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

% Repeat
%-----------------------------------------------------------------------
conf_class_repeat        = cfg_const;
conf_class_repeat.name   = 'Repeat';
conf_class_repeat.tag    = 'type';
conf_class_repeat.val    = {'cfg_repeat'};
conf_class_repeat.hidden = true;
conf_class_repeat.help   = {'Hidden field that gives the hint to cfg_struct2cfg which class to create.'};

%% Item generators

% Branch
%-----------------------------------------------------------------------
conf_branch      = cfg_exbranch;
conf_branch.name = 'Branch';
conf_branch.tag  = 'conf_branch';
conf_branch.val  = {conf_class_branch, conf_name, conf_tag, conf_val, conf_check, conf_help};
conf_branch.help = help2cell('cfg_branch');
conf_branch.prog = @cfg_cfg_pass;
conf_branch.vout = @cfg_cfg_vout;

% Choice
%-----------------------------------------------------------------------
conf_choice      = cfg_exbranch;
conf_choice.name = 'Choice';
conf_choice.tag  = 'conf_choice';
conf_choice.val  = {conf_class_choice, conf_name, conf_tag, conf_values, conf_check, conf_help};
conf_choice.help = help2cell('cfg_choice');
conf_choice.prog = @cfg_cfg_pass;
conf_choice.vout = @cfg_cfg_vout;

% Const
%-----------------------------------------------------------------------
conf_const      = cfg_exbranch;
conf_const.name = 'Const';
conf_const.tag  = 'conf_const';
conf_const.val  = {conf_class_const, conf_name, conf_tag, conf_val_single, ...
                   conf_check, conf_help, conf_def};
conf_const.help = help2cell('cfg_const');
conf_const.prog = @cfg_cfg_pass;
conf_const.vout = @cfg_cfg_vout;

% Entry
%-----------------------------------------------------------------------
conf_entry      = cfg_exbranch;
conf_entry.name = 'Entry';
conf_entry.tag  = 'conf_entry';
conf_entry.val  = {conf_class_entry, conf_name, conf_tag, conf_strtype, ...
                   conf_extras, conf_num_any, conf_check, conf_help, conf_def};
conf_entry.help = help2cell('cfg_entry');
conf_entry.prog = @cfg_cfg_pass;
conf_entry.vout = @cfg_cfg_vout;

% Exbranch
%-----------------------------------------------------------------------
conf_exbranch      = cfg_exbranch;
conf_exbranch.name = 'Exbranch';
conf_exbranch.tag  = 'conf_exbranch';
conf_exbranch.val  = {conf_class_exbranch, conf_name, conf_tag, conf_val, ...
                    conf_prog, conf_vout, conf_check, conf_help};
conf_exbranch.help = help2cell('cfg_exbranch');
conf_exbranch.prog = @cfg_cfg_pass;
conf_exbranch.vout = @cfg_cfg_vout;

% Files
%-----------------------------------------------------------------------
conf_files      = cfg_exbranch;
conf_files.name = 'Files';
conf_files.tag  = 'conf_files';
conf_files.val  = {conf_class_files, conf_name, conf_tag, conf_filter, ...
                   conf_ufilter, conf_dir, conf_num, conf_check, conf_help, conf_def};
conf_files.help = help2cell('cfg_files');
conf_files.prog = @cfg_cfg_pass;
conf_files.vout = @cfg_cfg_vout;

% Menu
%-----------------------------------------------------------------------
conf_menu      = cfg_exbranch;
conf_menu.name = 'Menu';
conf_menu.tag  = 'conf_menu';
conf_menu.val  = {conf_class_menu, conf_name, conf_tag, conf_labels, ...
                  conf_values, conf_check, conf_help, conf_def};
conf_menu.help = help2cell('cfg_menu');
conf_menu.prog = @cfg_cfg_pass;
conf_menu.vout = @cfg_cfg_vout;
conf_menu.check = @cfg_cfg_labels_values;

% repeat
%-----------------------------------------------------------------------
conf_repeat      = cfg_exbranch;
conf_repeat.name = 'Repeat';
conf_repeat.tag  = 'conf_repeat';
conf_repeat.val  = {conf_class_repeat, conf_name, conf_tag, conf_values, ...
                    conf_num, conf_forcestruct, conf_check, conf_help};
conf_repeat.help = help2cell('cfg_repeat');
conf_repeat.prog = @cfg_cfg_pass;
conf_repeat.vout = @cfg_cfg_vout;

%% Output nodes

% Generate code
%-----------------------------------------------------------------------

gencode_fname         = cfg_entry;
gencode_fname.name    = 'Output filename';
gencode_fname.tag     = 'gencode_fname';
gencode_fname.strtype = 's';
gencode_fname.num     = [1 Inf];
gencode_fname.help    = {'Filename for generated .m File.'};

gencode_dir         = cfg_files;
gencode_dir.name    = 'Output directory';
gencode_dir.tag     = 'gencode_dir';
gencode_dir.filter  = 'dir';
gencode_dir.num     = [1 1];
gencode_dir.help    = {'Output directory for generated .m File.'};

gencode_var         = cfg_entry;
gencode_var.name    = 'Root node of config';
gencode_var.tag     = 'gencode_var';
gencode_var.strtype = 'e';
gencode_var.num     = [1 1];
gencode_var.help    = {['This should be a dependency input from the root ' ...
                    'node of the application''s configuration tree. Use ' ...
                    'the output of the root configuration item directly, ' ...
                    'it is not necessary to run "Generate object tree" first.']};

gencode_o_def       = cfg_menu;
gencode_o_def.name  = 'Create Defaults File';
gencode_o_def.tag   = 'gencode_o_def';
gencode_o_def.labels= {'No',  'Yes'};
gencode_o_def.values= {false, true};
gencode_o_def.help  = {'The defaults file can be used to document all possible job structures. Its use to set all defaults is deprecated, .def function handles should be used instead.'};

gencode_o_mlb       = cfg_menu;
gencode_o_mlb.name  = 'Create mlbatch_appcfg File';
gencode_o_mlb.tag   = 'gencode_o_mlb';
gencode_o_mlb.labels= {'No',  'Yes'};
gencode_o_mlb.values= {false, true};
gencode_o_mlb.help  = {'The cfg_mlbatch_appcfg file can be used if the toolbox should be found by MATLABBATCH automatically.'};

gencode_o_path       = cfg_menu;
gencode_o_path.name  = 'Create Code for addpath()';
gencode_o_path.tag   = 'gencode_o_path';
gencode_o_path.labels= {'No',  'Yes'};
gencode_o_path.values= {false, true};
gencode_o_path.help  = {'If the toolbox resides in a non-MATLAB path, code can be generated to automatically add the configuration file to MATLAB path.'};

gencode_opts        = cfg_branch;
gencode_opts.name   = 'Options';
gencode_opts.tag    = 'gencode_opts';
gencode_opts.val    = {gencode_o_def gencode_o_mlb gencode_o_path};
gencode_opts.help   = {'Code generation options.'};

gencode_gen         = cfg_exbranch;
gencode_gen.name    = 'Generate code';
gencode_gen.tag     = 'gencode_gen';
gencode_gen.val     = {gencode_fname, gencode_dir, gencode_var, gencode_opts};
gencode_gen.help    = {['Generate code from a cfg_item tree. This tree can ' ...
                    'be either a struct (as returned from the ConfGUI ' ...
                    'modules) or a cfg_item object tree.']};
gencode_gen.prog    = @cfg_cfg_gencode;
gencode_gen.vout    = @vout_cfg_gencode;

% Generate object tree without code generation
%-----------------------------------------------------------------------

genobj_var         = cfg_entry;
genobj_var.name    = 'Root node of config';
genobj_var.tag     = 'genobj_var';
genobj_var.strtype = 'e';
genobj_var.num     = [1 1];
genobj_var.help    = {'This should be a dependency input from the root node of the application''s configuration tree.'};

genobj_gen         = cfg_exbranch;
genobj_gen.name    = 'Generate object tree';
genobj_gen.tag     = 'genobj_gen';
genobj_gen.val     = {genobj_var};
genobj_gen.help    = {['Generate a cfg_item tree as a variable. This can ' ...
                    'be useful to test a configuration before it is saved ' ...
                    'into files.']};
genobj_gen.prog    = @cfg_cfg_genobj;
genobj_gen.vout    = @vout_cfg_genobj;

%% Assemble Menu

% Data entry nodes
%-----------------------------------------------------------------------

menu_entry        = cfg_choice;
menu_entry.name   = 'Data entry items';
menu_entry.tag    = 'menu_entry';
menu_entry.values = {conf_entry, conf_files, conf_menu, conf_const};
menu_entry.help   = {'These items are used to enter data that will be passed to the computation code.'};

% Tree structuring nodes
%-----------------------------------------------------------------------

menu_struct        = cfg_choice;
menu_struct.name   = 'Tree structuring items';
menu_struct.tag    = 'menu_struct';
menu_struct.values = {conf_branch, conf_exbranch, conf_choice, conf_repeat};
menu_struct.help   = {'These items collect data entry items and build a menu structure.'};

% Root node
%-----------------------------------------------------------------------

menu_cfg        = cfg_choice;
menu_cfg.name   = 'ConfGUI';
menu_cfg.tag    = 'menu_cfg';
menu_cfg.values = {menu_entry, menu_struct, gencode_gen, genobj_gen};
menu_cfg.help   = help2cell(mfilename);

%% Helper functions

function out = cfg_cfg_genobj(varargin)
if isa(varargin{1}.genobj_var, 'cfg_item')
    % use object tree "as is"
    out.c0 = varargin{1}.gencode_var;
else
    % Transform struct into class based tree
    out.c0 = cfg_struct2cfg(varargin{1}.genobj_var);
end
[u1 out.djob]  = harvest(out.c0, out.c0, true, true);

function out = cfg_cfg_gencode(varargin)
if isa(varargin{1}.gencode_var, 'cfg_item')
    % use object tree "as is"
    out.c0 = varargin{1}.gencode_var;
else
    % Transform struct into class based tree
    out.c0 = cfg_struct2cfg(varargin{1}.gencode_var);
end
% Generate code
[str tag] = gencode(out.c0,'',{});
[p n e] = fileparts(varargin{1}.gencode_fname);
out.cfg_file{1} = fullfile(varargin{1}.gencode_dir{1}, [n '.m']);
fid = fopen(out.cfg_file{1}, 'wt');
fprintf(fid, 'function %s = %s\n', tag, n);
fprintf(fid, ...
        ['%% ''%s'' - MATLABBATCH configuration\n' ...
         '%% This MATLABBATCH configuration file has been generated automatically\n' ...
         '%% by MATLABBATCH using ConfGUI. It describes menu structure, validity\n' ...
         '%% constraints and links to run time code.\n' ...
         '%% Changes to this file will be overwritten if the ConfGUI batch is executed again.\n' ...
         '%% Created at %s.\n'], out.c0.name, datestr(now, 31));
fprintf(fid, '%s\n', str{:});
if varargin{1}.gencode_opts.gencode_o_path
    fprintf(fid, '%% ---------------------------------------------------------------------\n');
    fprintf(fid, '%% add path to this mfile\n');
    fprintf(fid, '%% ---------------------------------------------------------------------\n');
    fprintf(fid, 'addpath(fileparts(mfilename(''fullpath'')));\n');
end
fclose(fid);
if varargin{1}.gencode_opts.gencode_o_def
    % Generate defaults file
    [u1 out.djob]  = harvest(out.c0, out.c0, true, true);
    [str dtag] = gencode(out.djob, sprintf('%s_def', tag));
    dn = sprintf('%s_def', n);
    out.def_file{1} = fullfile(varargin{1}.gencode_dir{1}, sprintf('%s.m', dn));
    fid = fopen(out.def_file{1}, 'wt');
    fprintf(fid, 'function %s = %s\n', dtag, dn);
    fprintf(fid, ...
        ['%% ''%s'' - MATLABBATCH defaults\n' ...
        '%% This MATLABBATCH defaults file has been generated automatically\n' ...
        '%% by MATLABBATCH using ConfGUI. It contains all pre-defined values for\n' ...
        '%% menu items and provides a full documentation of all fields that may\n' ...
        '%% be present in a job variable for this application.\n' ...
        '%% Changes to this file will be overwritten if the ConfGUI batch is executed again.\n' ...
        '%% Created at %s.\n'], out.c0.name, datestr(now, 31));
    fprintf(fid, '%s\n', str{:});
    fclose(fid);
end
if varargin{1}.gencode_opts.gencode_o_mlb
    % Generate cfg_util initialisation file
    out.mlb_file{1} = fullfile(varargin{1}.gencode_dir{1}, 'cfg_mlbatch_appcfg.m');
    fid = fopen(out.mlb_file{1}, 'wt');
    fprintf(fid, 'function [cfg, def] = cfg_mlbatch_appcfg(varargin)\n');
    fprintf(fid, ...
        ['%% ''%s'' - MATLABBATCH cfg_util initialisation\n' ...
        '%% This MATLABBATCH initialisation file can be used to load application\n' ...
        '%%              ''%s''\n' ...
        '%% into cfg_util. This can be done manually by running this file from\n' ...
        '%% MATLAB command line or automatically when cfg_util is initialised.\n' ...
        '%% The directory containing this file and the configuration file\n' ...
        '%%              ''%s''\n' ...
        '%% must be in MATLAB''s path variable.\n' ...
        '%% Created at %s.\n\n'], ...
        out.c0.name, out.c0.name, n, datestr(now, 31));

    fprintf(fid, 'if ~isdeployed\n');
    fprintf(fid, '    %% Get path to this file and add it to MATLAB path.\n');
    fprintf(fid, ['    %% If the configuration file is stored in another place, the ' ...
        'path must be adjusted here.\n']);
    fprintf(fid, '    p = fileparts(mfilename(''fullpath''));\n');
    fprintf(fid, '    addpath(p);\n');
    fprintf(fid, 'end\n');
    fprintf(fid, '%% run configuration main & def function, return output\n');
    fprintf(fid, 'cfg = %s;\n', n);
    if varargin{1}.gencode_opts.gencode_o_def
        fprintf(fid, 'def = %s;\n', dn);
    else
        fprintf(fid, 'def = [];\n');
    end
    fclose(fid);
end

function out = cfg_cfg_pass(varargin)
% just pass input to output
out = varargin{1};

function str = cfg_cfg_labels_values(varargin)
% Check whether a menu has the same number of labels and values items
if numel(varargin{1}.labels) == numel(varargin{1}.values)
    str = '';
else
    str = 'Number of labels must match number of values.';
end

function vout = cfg_cfg_vout(varargin)
% cfg_struct2cfg returns its output immediately, so a subscript '(1)' is
% appropriate.
vout = cfg_dep;
vout.sname = sprintf('%s (%s)', varargin{1}.name, varargin{1}.type);
vout.src_output = substruct('()', {1});

function vout = vout_cfg_genobj(varargin)
vout(1)            = cfg_dep;
vout(1).sname      = 'Configuration Object Tree';
vout(1).src_output = substruct('.','c0');
vout(1).tgt_spec   = cfg_findspec({{'strtype','e'}});
vout(2)            = cfg_dep;
vout(2).sname      = 'Configuration Defaults Variable';
vout(2).src_output = substruct('.','djob');
vout(2).tgt_spec   = cfg_findspec({{'strtype','e'}});

function vout = vout_cfg_gencode(varargin)
vout(1)            = cfg_dep;
vout(1).sname      = 'Generated Configuration File';
vout(1).src_output = substruct('.', 'cfg_file');
vout(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
vout(2)            = cfg_dep;
vout(2).sname      = 'Configuration Object Tree';
vout(2).src_output = substruct('.','c0');
vout(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
if islogical(varargin{1}.gencode_opts.gencode_o_def) && varargin{1}.gencode_opts.gencode_o_def
    vout(3)            = cfg_dep;
    vout(3).sname      = 'Generated Defaults File';
    vout(3).src_output = substruct('.', 'def_file');
    vout(3).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
    vout(4)            = cfg_dep;
    vout(4).sname      = 'Configuration Defaults Variable';
    vout(4).src_output = substruct('.','djob');
    vout(4).tgt_spec   = cfg_findspec({{'strtype','e'}});
end
if islogical(varargin{1}.gencode_opts.gencode_o_mlb) && varargin{1}.gencode_opts.gencode_o_mlb
    vout(end+1)            = cfg_dep;
    vout(end).sname      = 'Generated Initialisation File';
    vout(end).src_output = substruct('.', 'mlb_file');
    vout(end).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
end
