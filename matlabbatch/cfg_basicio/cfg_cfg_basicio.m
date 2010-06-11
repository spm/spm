function cfg_basicio = cfg_cfg_basicio
% 'BasicIO' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2010-06-11 13:56:06.
% ---------------------------------------------------------------------
% files Files to move/copy/delete
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'Files to move/copy/delete';
files.help    = {'These files will be moved, copied or deleted.'};
files.filter = 'any';
files.ufilter = '.*';
files.num     = [0 Inf];
% ---------------------------------------------------------------------
% moveto Move to
% ---------------------------------------------------------------------
moveto         = cfg_files;
moveto.tag     = 'moveto';
moveto.name    = 'Move to';
moveto.help    = {'Files will be moved to the specified directory.'};
moveto.filter = 'dir';
moveto.ufilter = '.*';
moveto.num     = [1 1];
% ---------------------------------------------------------------------
% copyto Copy to
% ---------------------------------------------------------------------
copyto         = cfg_files;
copyto.tag     = 'copyto';
copyto.name    = 'Copy to';
copyto.help    = {'Files will be moved to the specified directory.'};
copyto.filter = 'dir';
copyto.ufilter = '.*';
copyto.num     = [1 1];
% ---------------------------------------------------------------------
% moveto Move to
% ---------------------------------------------------------------------
moveto1         = cfg_files;
moveto1.tag     = 'moveto';
moveto1.name    = 'Move to';
moveto1.help    = {'Files will be moved to the specified directory.'};
moveto1.filter = 'dir';
moveto1.ufilter = '.*';
moveto1.num     = [1 1];
% ---------------------------------------------------------------------
% pattern Pattern
% ---------------------------------------------------------------------
pattern         = cfg_entry;
pattern.tag     = 'pattern';
pattern.name    = 'Pattern';
pattern.help    = {'The regular expression pattern to look for.'};
pattern.strtype = 's';
pattern.num     = [1  Inf];
% ---------------------------------------------------------------------
% repl Replacement
% ---------------------------------------------------------------------
repl         = cfg_entry;
repl.tag     = 'repl';
repl.name    = 'Replacement';
repl.help    = {'This string (or pattern) will be inserted instead.'};
repl.strtype = 's';
repl.num     = [1  Inf];
% ---------------------------------------------------------------------
% patrep Pattern/Replacement Pair
% ---------------------------------------------------------------------
patrep         = cfg_branch;
patrep.tag     = 'patrep';
patrep.name    = 'Pattern/Replacement Pair';
patrep.val     = {pattern repl };
% ---------------------------------------------------------------------
% patreplist Pattern/Replacement List
% ---------------------------------------------------------------------
patreplist         = cfg_repeat;
patreplist.tag     = 'patreplist';
patreplist.name    = 'Pattern/Replacement List';
patreplist.help    = {'Regexprep supports a list of multiple patterns and corresponding replacements. These will be applied to the filename portion (without path, without extension) one after another. E.g., if your filename is ''testimage(.nii)'', and you replace ''test'' with ''xyz'' and ''xyzim'' with ''newtestim'', the final filename will be ''newtestimage.nii''.'};
patreplist.values  = {patrep };
patreplist.num     = [1 Inf];
% ---------------------------------------------------------------------
% unique Unique Filenames
% ---------------------------------------------------------------------
unique         = cfg_menu;
unique.tag     = 'unique';
unique.name    = 'Unique Filenames';
unique.help    = {
                  'If the regexprep operation results in identical output filenames for two or more input files, these can not be written/renamed to their new location without loosing data. If you are sure that your regexprep patterns produce unique filenames, you do not need to care about this.'
                  'If you choose to append a running number, it will be zero-padded to make sure alphabetical sort of filenames returns them in the same order as the input files are.'
                  }';
unique.labels = {
                 'Don''t Care'
                 'Append Index Number'
                 }';
unique.values = {
                 false
                 true
                 }';
% ---------------------------------------------------------------------
% moveren Move and Rename
% ---------------------------------------------------------------------
moveren         = cfg_branch;
moveren.tag     = 'moveren';
moveren.name    = 'Move and Rename';
moveren.val     = {moveto1 patreplist unique };
moveren.help    = {'The input files will be moved to the specified target folder. In addition, their filenames (not paths, not extensions) will be changed by replacing regular expression patterns using MATLABs regexprep function. Please consult MATLAB help and HTML documentation for how to specify regular expressions.'};
% ---------------------------------------------------------------------
% copyto Copy to
% ---------------------------------------------------------------------
copyto1         = cfg_files;
copyto1.tag     = 'copyto';
copyto1.name    = 'Copy to';
copyto1.help    = {'Files will be moved to the specified directory.'};
copyto1.filter = 'dir';
copyto1.ufilter = '.*';
copyto1.num     = [1 1];
% ---------------------------------------------------------------------
% pattern Pattern
% ---------------------------------------------------------------------
pattern         = cfg_entry;
pattern.tag     = 'pattern';
pattern.name    = 'Pattern';
pattern.help    = {'The regular expression pattern to look for.'};
pattern.strtype = 's';
pattern.num     = [1  Inf];
% ---------------------------------------------------------------------
% repl Replacement
% ---------------------------------------------------------------------
repl         = cfg_entry;
repl.tag     = 'repl';
repl.name    = 'Replacement';
repl.help    = {'This string (or pattern) will be inserted instead.'};
repl.strtype = 's';
repl.num     = [1  Inf];
% ---------------------------------------------------------------------
% patrep Pattern/Replacement Pair
% ---------------------------------------------------------------------
patrep         = cfg_branch;
patrep.tag     = 'patrep';
patrep.name    = 'Pattern/Replacement Pair';
patrep.val     = {pattern repl };
% ---------------------------------------------------------------------
% patreplist Pattern/Replacement List
% ---------------------------------------------------------------------
patreplist         = cfg_repeat;
patreplist.tag     = 'patreplist';
patreplist.name    = 'Pattern/Replacement List';
patreplist.help    = {'Regexprep supports a list of multiple patterns and corresponding replacements. These will be applied to the filename portion (without path, without extension) one after another. E.g., if your filename is ''testimage(.nii)'', and you replace ''test'' with ''xyz'' and ''xyzim'' with ''newtestim'', the final filename will be ''newtestimage.nii''.'};
patreplist.values  = {patrep };
patreplist.num     = [1 Inf];
% ---------------------------------------------------------------------
% unique Unique Filenames
% ---------------------------------------------------------------------
unique         = cfg_menu;
unique.tag     = 'unique';
unique.name    = 'Unique Filenames';
unique.help    = {
                  'If the regexprep operation results in identical output filenames for two or more input files, these can not be written/renamed to their new location without loosing data. If you are sure that your regexprep patterns produce unique filenames, you do not need to care about this.'
                  'If you choose to append a running number, it will be zero-padded to make sure alphabetical sort of filenames returns them in the same order as the input files are.'
                  }';
unique.labels = {
                 'Don''t Care'
                 'Append Index Number'
                 }';
unique.values = {
                 false
                 true
                 }';
% ---------------------------------------------------------------------
% copyren Copy and Rename
% ---------------------------------------------------------------------
copyren         = cfg_branch;
copyren.tag     = 'copyren';
copyren.name    = 'Copy and Rename';
copyren.val     = {copyto1 patreplist unique };
copyren.help    = {'The input files will be copied to the specified target folder. In addition, their filenames (not paths, not extensions) will be changed by replacing regular expression patterns using MATLABs regexprep function. Please consult MATLAB help and HTML documentation for how to specify regular expressions.'};
% ---------------------------------------------------------------------
% delete Delete
% ---------------------------------------------------------------------
delete         = cfg_const;
delete.tag     = 'delete';
delete.name    = 'Delete';
delete.val = {false};
delete.help    = {'The selected files will be deleted.'};
% ---------------------------------------------------------------------
% action Action
% ---------------------------------------------------------------------
action         = cfg_choice;
action.tag     = 'action';
action.name    = 'Action';
action.values  = {moveto copyto moveren copyren delete };
% ---------------------------------------------------------------------
% file_move Move/Delete Files
% ---------------------------------------------------------------------
file_move         = cfg_exbranch;
file_move.tag     = 'file_move';
file_move.name    = 'Move/Delete Files';
file_move.val     = {files action };
file_move.help    = {'Move or delete files.'};
file_move.prog = @cfg_run_file_move;
file_move.vout = @cfg_vout_file_move;
% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'New working directory.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];
% ---------------------------------------------------------------------
% cfg_cd Change Directory
% ---------------------------------------------------------------------
cfg_cd         = cfg_exbranch;
cfg_cd.tag     = 'cfg_cd';
cfg_cd.name    = 'Change Directory';
cfg_cd.val     = {dir };
cfg_cd.help    = {'Change working directory.'};
cfg_cd.prog = @cfg_run_cd;
% ---------------------------------------------------------------------
% parent Parent Directory
% ---------------------------------------------------------------------
parent         = cfg_files;
parent.tag     = 'parent';
parent.name    = 'Parent Directory';
parent.help    = {'Directory where the new directory will be created.'};
parent.filter = 'dir';
parent.ufilter = '.*';
parent.num     = [1 1];
% ---------------------------------------------------------------------
% name New Directory Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'New Directory Name';
name.help    = {'Name for the new directory.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% cfg_mkdir Make Directory
% ---------------------------------------------------------------------
cfg_mkdir         = cfg_exbranch;
cfg_mkdir.tag     = 'cfg_mkdir';
cfg_mkdir.name    = 'Make Directory';
cfg_mkdir.val     = {parent name };
cfg_mkdir.help    = {'Create a new directory.'};
cfg_mkdir.prog = @cfg_run_mkdir;
cfg_mkdir.vout = @cfg_vout_mkdir;
% ---------------------------------------------------------------------
% name Input Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Input Name';
name.help    = {'Enter a name for this directory selection. This name will be displayed in the ''Dependency'' listing as output name.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% dirs Directory
% ---------------------------------------------------------------------
dirs1         = cfg_files;
dirs1.tag     = 'dirs';
dirs1.name    = 'Directory';
dirs1.help    = {'Select a directory.'};
dirs1.filter = 'dir';
dirs1.ufilter = '.*';
dirs1.num     = [1 1];
% ---------------------------------------------------------------------
% dirs Directories
% ---------------------------------------------------------------------
dirs         = cfg_repeat;
dirs.tag     = 'dirs';
dirs.name    = 'Directories';
dirs.help    = {'Select one or more directories.'};
dirs.values  = {dirs1 };
dirs.num     = [1 Inf];
% ---------------------------------------------------------------------
% cfg_named_dir Named Directory Selector
% ---------------------------------------------------------------------
cfg_named_dir         = cfg_exbranch;
cfg_named_dir.tag     = 'cfg_named_dir';
cfg_named_dir.name    = 'Named Directory Selector';
cfg_named_dir.val     = {name dirs };
cfg_named_dir.help    = {'Named Directory Selector allows to select directories that can be referenced as common input by other modules.'};
cfg_named_dir.prog = @cfg_run_named_dir;
cfg_named_dir.vout = @cfg_vout_named_dir;
% ---------------------------------------------------------------------
% name Input Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Input Name';
name.help    = {'Enter a name for this file selection. This name will be displayed in the ''Dependency'' listing as output name.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% files File Set
% ---------------------------------------------------------------------
files1         = cfg_files;
files1.tag     = 'files';
files1.name    = 'File Set';
files1.help    = {'Select a set of files.'};
files1.filter = 'any';
files1.ufilter = '.*';
files1.num     = [0 Inf];
% ---------------------------------------------------------------------
% files File Sets
% ---------------------------------------------------------------------
files         = cfg_repeat;
files.tag     = 'files';
files.name    = 'File Sets';
files.help    = {'Select one or more sets of files. Each set can be passed separately as a dependency to other modules.'};
files.values  = {files1 };
files.num     = [1 Inf];
% ---------------------------------------------------------------------
% cfg_named_file Named File Selector
% ---------------------------------------------------------------------
cfg_named_file         = cfg_exbranch;
cfg_named_file.tag     = 'cfg_named_file';
cfg_named_file.name    = 'Named File Selector';
cfg_named_file.val     = {name files };
cfg_named_file.help    = {
                          'Named File Selector allows to select sets of files that can be referenced as common input by other modules.'
                          'In addition to file outputs, an index vector is provided that can be used to separate the files again using ''File Set Split'' module. This can be useful if the same sets of files have to be processed separately in one module and as a single set in another.'
                          }';
cfg_named_file.prog = @cfg_run_named_file;
cfg_named_file.vout = @cfg_vout_named_file;
% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Directory to select files from.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 Inf];
% ---------------------------------------------------------------------
% filter Filter
% ---------------------------------------------------------------------
filter         = cfg_entry;
filter.tag     = 'filter';
filter.name    = 'Filter';
filter.help    = {'A regular expression to filter files.'};
filter.strtype = 's';
filter.num     = [1  Inf];
% ---------------------------------------------------------------------
% rec Descend into subdirectories
% ---------------------------------------------------------------------
rec         = cfg_menu;
rec.tag     = 'rec';
rec.name    = 'Descend into subdirectories';
rec.help    = {'Files can be selected from the specified directory only or from the specified directory and all its subdirectories.'};
rec.labels = {
              'Yes'
              'No'
              }';
rec.values = {
              'FPListRec'
              'FPList'
              }';
% ---------------------------------------------------------------------
% file_fplist File Selector (Batch Mode)
% ---------------------------------------------------------------------
file_fplist         = cfg_exbranch;
file_fplist.tag     = 'file_fplist';
file_fplist.name    = 'File Selector (Batch Mode)';
file_fplist.val     = {dir filter rec };
file_fplist.help    = {'Select files from a directory using cfg_getfile(''FPList'',...).'};
file_fplist.prog = @cfg_run_file_fplist;
file_fplist.vout = @cfg_vout_file_fplist;
% ---------------------------------------------------------------------
% files Files
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'Files';
files.help    = {'Files to be filtered.'};
files.filter = 'any';
files.ufilter = '.*';
files.num     = [1 Inf];
% ---------------------------------------------------------------------
% typ Typ
% ---------------------------------------------------------------------
typ         = cfg_entry;
typ.tag     = 'typ';
typ.name    = 'Typ';
typ.help    = {'Allowed types are (see cfg_getfile): ''any'', ''image'', ''xml'', ''mat'', ''batch'', ''dir'' or a regular expression.'};
typ.strtype = 's';
typ.num     = [1  Inf];
% ---------------------------------------------------------------------
% filter Filter
% ---------------------------------------------------------------------
filter         = cfg_entry;
filter.tag     = 'filter';
filter.name    = 'Filter';
filter.help    = {'A regular expression to filter files (applied after filtering for ''Typ'').'};
filter.strtype = 's';
filter.num     = [1  Inf];
% ---------------------------------------------------------------------
% frames Frames
% ---------------------------------------------------------------------
frames         = cfg_entry;
frames.tag     = 'frames';
frames.name    = 'Frames';
frames.help    = {'This option is currently completely ignored!'};
frames.strtype = 's';
frames.num     = [0  Inf];
% ---------------------------------------------------------------------
% file_filter File Filter
% ---------------------------------------------------------------------
file_filter         = cfg_exbranch;
file_filter.tag     = 'file_filter';
file_filter.name    = 'File Filter';
file_filter.val     = {files typ filter frames };
file_filter.help    = {'Filter a list of files using cfg_getfile.'};
file_filter.prog = @cfg_run_file_filter;
file_filter.vout = @cfg_vout_file_filter;
% ---------------------------------------------------------------------
% name File Set Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'File Set Name';
name.help    = {'Enter a name for this file selection. This name will be displayed in the ''Dependency'' listing as output name.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% files Input File Set
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'Input File Set';
files.help    = {'The input file set will be split at the indices given in the ''#files per set'' collection.'};
files.filter = 'any';
files.ufilter = '.*';
files.num     = [1 Inf];
% ---------------------------------------------------------------------
% index Selection Index
% ---------------------------------------------------------------------
index         = cfg_entry;
index.tag     = 'index';
index.name    = 'Selection Index';
index.help    = {
                 'Enter the index vector of files belonging to this set. E.g. to select files 1 to 10 from the input file set, enter [1:10].'
                 'To split the combined list of all files selected by a "Named File Selector", enter the index vectors here as dependency.'
                 }';
index.strtype = 'n';
index.num     = [1  Inf];
% ---------------------------------------------------------------------
% filesets Output File Sets
% ---------------------------------------------------------------------
filesets         = cfg_repeat;
filesets.tag     = 'filesets';
filesets.name    = 'Output File Sets';
filesets.help    = {'For each file set to be created, enter an index vector that selects files from the input list. An additional output file set will be created matching files not selected by any index.'};
filesets.values  = {index };
filesets.num     = [1 Inf];
% ---------------------------------------------------------------------
% cfg_file_split File Set Split
% ---------------------------------------------------------------------
cfg_file_split         = cfg_exbranch;
cfg_file_split.tag     = 'cfg_file_split';
cfg_file_split.name    = 'File Set Split';
cfg_file_split.val     = {name files filesets };
cfg_file_split.help    = {
                          'Split a list of files into multiple parts.'
                          'With this utility, a list of files can be split into parts based on index vectors. Each index vector is a list of numbers where number N selects the Nth file from the input list. Index vectors may overlap, so files may belong to more than one part. Files not matching any index will be returned in a separate list. If an index is larger than the number of files passed it will be ignored.'
                          }';
cfg_file_split.prog = @cfg_run_file_split;
cfg_file_split.vout = @cfg_vout_file_split;
% ---------------------------------------------------------------------
% name Input Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Input Name';
name.help    = {'Enter a name for this variable. This name will be displayed in the ''Dependency'' listing as output name.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% input Input Variable
% ---------------------------------------------------------------------
input         = cfg_entry;
input.tag     = 'input';
input.name    = 'Input Variable';
input.help    = {'Enter a MATLAB variable. This can be a variable in the MATLAB workspace, or any other valid MATLAB statement which evaluates to a single variable.'};
input.strtype = 'e';
input.num     = [];
% ---------------------------------------------------------------------
% cfg_named_input Named Input
% ---------------------------------------------------------------------
cfg_named_input         = cfg_exbranch;
cfg_named_input.tag     = 'cfg_named_input';
cfg_named_input.name    = 'Named Input';
cfg_named_input.val     = {name input };
cfg_named_input.help    = {'Named Input allows to enter any kind of MATLAB variable. This variable can be referenced as common input by other modules.'};
cfg_named_input.prog = @cfg_run_named_input;
cfg_named_input.vout = @cfg_vout_named_input;
% ---------------------------------------------------------------------
% matname .mat Filename
% ---------------------------------------------------------------------
matname         = cfg_files;
matname.tag     = 'matname';
matname.name    = '.mat Filename';
matname.help    = {'The name of the .mat file to load.'};
matname.filter = 'mat';
matname.ufilter = '.*';
matname.num     = [1 1];
% ---------------------------------------------------------------------
% allvars All Variables
% ---------------------------------------------------------------------
allvars         = cfg_const;
allvars.tag     = 'allvars';
allvars.name    = 'All Variables';
allvars.val = {true};
allvars.help    = {'Load all variables found in the .mat file.'};
% ---------------------------------------------------------------------
% varname Variable Name
% ---------------------------------------------------------------------
varname         = cfg_entry;
varname.tag     = 'varname';
varname.name    = 'Variable Name';
varname.check   = @(job)cfg_load_vars('check','isvarname',job);
varname.help    = {'Enter the name of the variable to be loaded from the .mat file.'};
varname.strtype = 's';
varname.num     = [1  Inf];
% ---------------------------------------------------------------------
% varnames Specified Variables
% ---------------------------------------------------------------------
varnames         = cfg_repeat;
varnames.tag     = 'varnames';
varnames.name    = 'Specified Variables';
varnames.help    = {'Enter a list of variable names to be loaded from the .mat file.'};
varnames.values  = {varname };
varnames.num     = [1 Inf];
% ---------------------------------------------------------------------
% loadvars Variables to load
% ---------------------------------------------------------------------
loadvars         = cfg_choice;
loadvars.tag     = 'loadvars';
loadvars.name    = 'Variables to load';
loadvars.help    = {'Choose whether all variables or a list of variables with specified names should be loaded.'};
loadvars.values  = {allvars varnames };
% ---------------------------------------------------------------------
% load_vars Load Variables from .mat File
% ---------------------------------------------------------------------
load_vars         = cfg_exbranch;
load_vars.tag     = 'load_vars';
load_vars.name    = 'Load Variables from .mat File';
load_vars.val     = {matname loadvars };
load_vars.help    = {'This function loads variables from a .mat file and passes them on as a dependency. It can load either all variables from a .mat file or a list of specified variables. In the first case, it will return a single struct variable. The variable names in the .mat file will become field names in this struct. If a list of variable names is given, each of them will be loaded into a separate variable.'};
load_vars.prog = @(job)cfg_load_vars('run',job);
load_vars.vout = @(job)cfg_load_vars('vout',job);
% ---------------------------------------------------------------------
% name Output Filename
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Output Filename';
name.help    = {'Output filename without any directory names.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% outdir Output Directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.help    = {'Directory where the file will be saved. Any directory components in the output filename will be stripped off and only this directory determines the path to the file.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% vname Variable Name
% ---------------------------------------------------------------------
vname         = cfg_entry;
vname.tag     = 'vname';
vname.name    = 'Variable Name';
vname.check   = @cfg_check_assignin;
vname.help    = {'Name for the variable in output file/struct. This must be a valid MATLAB variable name.'};
vname.strtype = 's';
vname.num     = [1  Inf];
% ---------------------------------------------------------------------
% vcont Variable Contents
% ---------------------------------------------------------------------
vcont         = cfg_entry;
vcont.tag     = 'vcont';
vcont.name    = 'Variable Contents';
vcont.help    = {'Contents to be saved. This can be any MATLAB variable or an output from another module, passed as dependency.'};
vcont.strtype = 'e';
vcont.num     = [];
% ---------------------------------------------------------------------
% vars Variable
% ---------------------------------------------------------------------
vars1         = cfg_branch;
vars1.tag     = 'vars';
vars1.name    = 'Variable';
vars1.val     = {vname vcont };
vars1.help    = {'For each variable, a name and its contents has to be specified.'};
% ---------------------------------------------------------------------
% vars Variables
% ---------------------------------------------------------------------
vars         = cfg_repeat;
vars.tag     = 'vars';
vars.name    = 'Variables';
vars.help    = {'Any number of variables can be saved in this file together.'};
vars.values  = {vars1 };
vars.num     = [1 Inf];
% ---------------------------------------------------------------------
% saveasstruct Save as
% ---------------------------------------------------------------------
saveasstruct         = cfg_menu;
saveasstruct.tag     = 'saveasstruct';
saveasstruct.name    = 'Save as';
saveasstruct.help    = {'Variables can be saved into the file individually or as a single struct variable, with the names of the variables used as fieldnames.'};
saveasstruct.labels = {
                       'Individual Variables'
                       'Struct Variable'
                       }';
saveasstruct.values = {
                       false
                       true
                       }';
% ---------------------------------------------------------------------
% cfg_save_vars Save Variables
% ---------------------------------------------------------------------
cfg_save_vars         = cfg_exbranch;
cfg_save_vars.tag     = 'cfg_save_vars';
cfg_save_vars.name    = 'Save Variables';
cfg_save_vars.val     = {name outdir vars saveasstruct };
cfg_save_vars.help    = {'Save a collection of variables to a .mat file.'};
cfg_save_vars.prog = @cfg_run_save_vars;
cfg_save_vars.vout = @cfg_vout_save_vars;
% ---------------------------------------------------------------------
% input Input variable
% ---------------------------------------------------------------------
input         = cfg_entry;
input.tag     = 'input';
input.name    = 'Input variable';
input.help    = {'Enter the input variable.'};
input.strtype = 'e';
input.num     = [Inf  Inf];
% ---------------------------------------------------------------------
% subsfield Field reference
% ---------------------------------------------------------------------
subsfield         = cfg_entry;
subsfield.tag     = 'subsfield';
subsfield.name    = 'Field reference';
subsfield.strtype = 's';
subsfield.num     = [1  Inf];
% ---------------------------------------------------------------------
% subsindc Cell index
% ---------------------------------------------------------------------
subsindc1         = cfg_entry;
subsindc1.tag     = 'subsindc';
subsindc1.name    = 'Cell index';
subsindc1.check   = @(job)cfg_run_subsrefvar('check','subsind',job);
subsindc1.help    = {'Enter a range of positive numbers or the string '':'' to reference all items.'};
subsindc1.strtype = 'e';
subsindc1.num     = [1  Inf];
% ---------------------------------------------------------------------
% subsindc Cell array reference
% ---------------------------------------------------------------------
subsindc         = cfg_repeat;
subsindc.tag     = 'subsindc';
subsindc.name    = 'Cell array reference';
subsindc.help    = {'For each dimension of an array, enter the subscripts to be indexed. Any multidimensional array can also be indexed by a single linear index.'};
subsindc.values  = {subsindc1 };
subsindc.num     = [1 Inf];
% ---------------------------------------------------------------------
% subsinda Array index
% ---------------------------------------------------------------------
subsinda1         = cfg_entry;
subsinda1.tag     = 'subsinda';
subsinda1.name    = 'Array index';
subsinda1.check   = @(job)cfg_run_subsrefvar('check','subsind',job);
subsinda1.help    = {'Enter a range of positive numbers or the string '':'' to reference all items.'};
subsinda1.strtype = 'e';
subsinda1.num     = [1  Inf];
% ---------------------------------------------------------------------
% subsinda Array reference
% ---------------------------------------------------------------------
subsinda         = cfg_repeat;
subsinda.tag     = 'subsinda';
subsinda.name    = 'Array reference';
subsinda.help    = {'For each dimension of an array, enter the subscripts to be indexed. Any multidimensional array can also be indexed by a single linear index.'};
subsinda.values  = {subsinda1 };
subsinda.num     = [1 Inf];
% ---------------------------------------------------------------------
% subsreference Part of variable to access
% ---------------------------------------------------------------------
subsreference         = cfg_repeat;
subsreference.tag     = 'subsreference';
subsreference.name    = 'Part of variable to access';
subsreference.help    = {
                         'Enter the sequence of subscript references. The leftmost reference is on the top of the list, the rightmost at the bottom.'
                         'To reference e.g. the field ''xY'' in a SPM variable, enter ''xY'' as ''Field reference''. To reference SPM.xY(1), enter ''xY'' as field reference, followed by 1 as an ''Array reference''.'
                         }';
subsreference.values  = {subsfield subsindc subsinda };
subsreference.num     = [1 Inf];
subsreference.forcestruct = true;
% ---------------------------------------------------------------------
% s String
% ---------------------------------------------------------------------
s         = cfg_const;
s.tag     = 's';
s.name    = 'String';
s.val{1}{1}.name = 'strtype';
s.val{1}{1}.value = 's';
% ---------------------------------------------------------------------
% n Natural Number
% ---------------------------------------------------------------------
n         = cfg_const;
n.tag     = 'n';
n.name    = 'Natural Number';
n.val{1}{1}.name = 'strtype';
n.val{1}{1}.value = 'n';
% ---------------------------------------------------------------------
% w Whole Number
% ---------------------------------------------------------------------
w         = cfg_const;
w.tag     = 'w';
w.name    = 'Whole Number';
w.val{1}{1}.name = 'strtype';
w.val{1}{1}.value = 'w';
% ---------------------------------------------------------------------
% i Integer
% ---------------------------------------------------------------------
i         = cfg_const;
i.tag     = 'i';
i.name    = 'Integer';
i.val{1}{1}.name = 'strtype';
i.val{1}{1}.value = 'i';
% ---------------------------------------------------------------------
% r Real Number
% ---------------------------------------------------------------------
r         = cfg_const;
r.tag     = 'r';
r.name    = 'Real Number';
r.val{1}{1}.name = 'strtype';
r.val{1}{1}.value = 'r';
% ---------------------------------------------------------------------
% f Function Handle
% ---------------------------------------------------------------------
f         = cfg_const;
f.tag     = 'f';
f.name    = 'Function Handle';
f.val{1}{1}.name = 'strtype';
f.val{1}{1}.value = 'f';
% ---------------------------------------------------------------------
% e Other Variable
% ---------------------------------------------------------------------
e         = cfg_const;
e.tag     = 'e';
e.name    = 'Other Variable';
e.val{1}{1}.name = 'strtype';
e.val{1}{1}.value = 'e';
% ---------------------------------------------------------------------
% nifti NIfTI Image File(s)
% ---------------------------------------------------------------------
nifti         = cfg_const;
nifti.tag     = 'nifti';
nifti.name    = 'NIfTI Image File(s)';
nifti.val{1}{1}.name = 'filter';
nifti.val{1}{1}.value = 'image';
% ---------------------------------------------------------------------
% mat .mat File(s)
% ---------------------------------------------------------------------
mat         = cfg_const;
mat.tag     = 'mat';
mat.name    = '.mat File(s)';
mat.val{1}{1}.name = 'filter';
mat.val{1}{1}.value = 'mat';
% ---------------------------------------------------------------------
% any Any File(s)
% ---------------------------------------------------------------------
any         = cfg_const;
any.tag     = 'any';
any.name    = 'Any File(s)';
any.val{1}{1}.name = 'filter';
any.val{1}{1}.value = 'any';
% ---------------------------------------------------------------------
% dir Directories
% ---------------------------------------------------------------------
dir         = cfg_const;
dir.tag     = 'dir';
dir.name    = 'Directories';
dir.val{1}{1}.name = 'filter';
dir.val{1}{1}.value = 'dir';
% ---------------------------------------------------------------------
% tgt_spec Type of referenced Variable
% ---------------------------------------------------------------------
tgt_spec         = cfg_choice;
tgt_spec.tag     = 'tgt_spec';
tgt_spec.name    = 'Type of referenced Variable';
tgt_spec.help    = {'This setting determines where the contents of the variable will be available as a dependency.'};
tgt_spec.values  = {s n w i r f e nifti mat any dir };
% ---------------------------------------------------------------------
% subsrefvar Access part of MATLAB variable
% ---------------------------------------------------------------------
subsrefvar         = cfg_exbranch;
subsrefvar.tag     = 'subsrefvar';
subsrefvar.name    = 'Access part of MATLAB variable';
subsrefvar.val     = {input subsreference tgt_spec };
subsrefvar.prog = @(job)cfg_run_subsrefvar('run',job);
subsrefvar.vout = @(job)cfg_run_subsrefvar('vout',job);
% ---------------------------------------------------------------------
% name Output Variable Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Output Variable Name';
name.check   = @cfg_check_assignin;
name.help    = {'Enter a valid MATLAB variable name. The contents of the input to "Output Item" will be assigned to a variable of this name in MATLAB workspace.'};
name.strtype = 's';
name.num     = [1  Inf];
% ---------------------------------------------------------------------
% output Output Item
% ---------------------------------------------------------------------
output         = cfg_entry;
output.tag     = 'output';
output.name    = 'Output Item';
output.help    = {'The contents that is passed to this input that will be assigned to the workspace variable whose name is given in "Output Variable Name".'};
output.strtype = 'e';
output.num     = [];
% ---------------------------------------------------------------------
% cfg_assignin Pass Output to Workspace
% ---------------------------------------------------------------------
cfg_assignin         = cfg_exbranch;
cfg_assignin.tag     = 'cfg_assignin';
cfg_assignin.name    = 'Pass Output to Workspace';
cfg_assignin.val     = {name output };
cfg_assignin.help    = {
                        'Assign a computation result to a workspace variable.'
                        'The value entered into "Output Item" will be assigned to a MATLAB workspace variable whose name is specified in "Output Variable Name". If this variable already exists, a new variable name will be generated. This can be useful to assess the results of computations with other MATLAB routines or for debugging.'
                        }';
cfg_assignin.prog = @cfg_run_assignin;
% ---------------------------------------------------------------------
% jobs Job File(s)
% ---------------------------------------------------------------------
jobs         = cfg_files;
jobs.tag     = 'jobs';
jobs.name    = 'Job File(s)';
jobs.help    = {'Select the job template(s). If multiple files are selected, they will be concatenated in selection order to form one job.'};
jobs.filter = 'batch';
jobs.ufilter = '.*';
jobs.num     = [1 Inf];
% ---------------------------------------------------------------------
% instr String
% ---------------------------------------------------------------------
instr         = cfg_entry;
instr.tag     = 'instr';
instr.name    = 'String';
instr.help    = {'Enter a string.'};
instr.strtype = 's';
instr.num     = [0  Inf];
% ---------------------------------------------------------------------
% ineval Evaluated Input
% ---------------------------------------------------------------------
ineval         = cfg_entry;
ineval.tag     = 'ineval';
ineval.name    = 'Evaluated Input';
ineval.help    = {'Enter an evaluated input.'};
ineval.strtype = 'e';
ineval.num     = [];
% ---------------------------------------------------------------------
% innifti NIfTI Images
% ---------------------------------------------------------------------
innifti         = cfg_files;
innifti.tag     = 'innifti';
innifti.name    = 'NIfTI Images';
innifti.help    = {'Select NIfTI Images'};
innifti.filter = 'image';
innifti.ufilter = '.*';
innifti.num     = [0 Inf];
% ---------------------------------------------------------------------
% inmat MATLAB .mat Files
% ---------------------------------------------------------------------
inmat         = cfg_files;
inmat.tag     = 'inmat';
inmat.name    = 'MATLAB .mat Files';
inmat.help    = {'Select MATLAB .mat files.'};
inmat.filter = 'mat';
inmat.ufilter = '.*';
inmat.num     = [0 Inf];
% ---------------------------------------------------------------------
% inany Any Files
% ---------------------------------------------------------------------
inany         = cfg_files;
inany.tag     = 'inany';
inany.name    = 'Any Files';
inany.help    = {'Select any kind of files.'};
inany.filter = 'any';
inany.ufilter = '.*';
inany.num     = [0 Inf];
% ---------------------------------------------------------------------
% indir Directory
% ---------------------------------------------------------------------
indir         = cfg_files;
indir.tag     = 'indir';
indir.name    = 'Directory';
indir.help    = {'Directory'};
indir.filter = 'dir';
indir.ufilter = '.*';
indir.num     = [1 Inf];
% ---------------------------------------------------------------------
% inputs Job Inputs
% ---------------------------------------------------------------------
inputs         = cfg_repeat;
inputs.tag     = 'inputs';
inputs.name    = 'Job Inputs';
inputs.help    = {'Assemble the set of input items for one run of the job.'};
inputs.values  = {instr ineval innifti inmat inany indir };
inputs.num     = [0 Inf];
% ---------------------------------------------------------------------
% runs Runs
% ---------------------------------------------------------------------
runs         = cfg_repeat;
runs.tag     = 'runs';
runs.name    = 'Runs';
runs.help    = {'Repeat "Job Inputs" for each run of the job, even if you want to specify no inputs in this batch itself. The number of "Job Inputs" items determines how often this job is run.'};
runs.values  = {inputs };
runs.num     = [1 Inf];
% ---------------------------------------------------------------------
% outstub Batch Filename Stub
% ---------------------------------------------------------------------
outstub         = cfg_entry;
outstub.tag     = 'outstub';
outstub.name    = 'Batch Filename Stub';
outstub.help    = {'The output filename will be generated by appending a job counter to this string.'};
outstub.strtype = 's';
outstub.num     = [1  Inf];
% ---------------------------------------------------------------------
% outdir Batch Directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Batch Directory';
outdir.help    = {'The generated batches will be saved into this folder.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% savejobs Save
% ---------------------------------------------------------------------
savejobs         = cfg_branch;
savejobs.tag     = 'savejobs';
savejobs.name    = 'Save';
savejobs.val     = {outstub outdir };
savejobs.help    = {'Specify filename stub and output directory to save the generated files.'};
% ---------------------------------------------------------------------
% dontsave Don't Save
% ---------------------------------------------------------------------
dontsave         = cfg_const;
dontsave.tag     = 'dontsave';
dontsave.name    = 'Don''t Save';
dontsave.val = {false};
dontsave.help    = {'Do not save the generated jobs.'};
% ---------------------------------------------------------------------
% save Save Generated Batch Jobs
% ---------------------------------------------------------------------
save         = cfg_choice;
save.tag     = 'save';
save.name    = 'Save Generated Batch Jobs';
save.help    = {'The generated batch jobs can be saved for reference or debugging purposes.'};
save.values  = {savejobs dontsave };
% ---------------------------------------------------------------------
% missing Missing Inputs
% ---------------------------------------------------------------------
missing         = cfg_menu;
missing.tag     = 'missing';
missing.name    = 'Missing Inputs';
missing.help    = {'Jobs with missing inputs (e.g. because of wrong input contents) can be skipped, while filled jobs are run. Alternatively, if any job has missing inputs, no jobs are run.'};
missing.labels = {
                  'Skip jobs with missing inputs, run filled jobs'
                  'Don''t run any jobs if missing inputs'
                  }';
missing.values = {
                  'skip'
                  'error'
                  }';
% ---------------------------------------------------------------------
% runjobs Run Batch Jobs
% ---------------------------------------------------------------------
runjobs         = cfg_exbranch;
runjobs.tag     = 'runjobs';
runjobs.name    = 'Run Batch Jobs';
runjobs.val     = {jobs runs save missing };
runjobs.help    = {'Load a set of job files, fill missing inputs and run the filled job. This automates the creation and execution of batch jobs for a large number of identical computations.'};
runjobs.prog = @cfg_run_runjobs;
runjobs.vout = @cfg_vout_runjobs;
% ---------------------------------------------------------------------
% evaluated Evaluated Input
% ---------------------------------------------------------------------
evaluated         = cfg_entry;
evaluated.tag     = 'evaluated';
evaluated.name    = 'Evaluated Input';
evaluated.strtype = 'e';
evaluated.num     = [];
% ---------------------------------------------------------------------
% string String
% ---------------------------------------------------------------------
string         = cfg_entry;
string.tag     = 'string';
string.name    = 'String';
string.strtype = 's';
string.num     = [];
% ---------------------------------------------------------------------
% anyfile Any File
% ---------------------------------------------------------------------
anyfile         = cfg_files;
anyfile.tag     = 'anyfile';
anyfile.name    = 'Any File';
anyfile.filter = 'any';
anyfile.ufilter = '.*';
anyfile.num     = [0 Inf];
% ---------------------------------------------------------------------
% images NIfTI Image(s)
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'NIfTI Image(s)';
images.filter = 'image';
images.ufilter = '.*';
images.num     = [0 Inf];
% ---------------------------------------------------------------------
% directory Directory
% ---------------------------------------------------------------------
directory         = cfg_files;
directory.tag     = 'directory';
directory.name    = 'Directory';
directory.filter = 'dir';
directory.ufilter = '.*';
directory.num     = [0 Inf];
% ---------------------------------------------------------------------
% inputs Inputs
% ---------------------------------------------------------------------
inputs         = cfg_repeat;
inputs.tag     = 'inputs';
inputs.name    = 'Inputs';
inputs.help    = {'Assemble the inputs to the called function in their correct order.'};
inputs.values  = {evaluated string anyfile images directory };
inputs.num     = [0 Inf];
inputs.forcestruct = true;
% ---------------------------------------------------------------------
% s String
% ---------------------------------------------------------------------
s         = cfg_const;
s.tag     = 's';
s.name    = 'String';
s.val = {true};
% ---------------------------------------------------------------------
% e Evaluated
% ---------------------------------------------------------------------
e         = cfg_const;
e.tag     = 'e';
e.name    = 'Evaluated';
e.val = {true};
% ---------------------------------------------------------------------
% n Natural number
% ---------------------------------------------------------------------
n         = cfg_const;
n.tag     = 'n';
n.name    = 'Natural number';
n.val = {true};
% ---------------------------------------------------------------------
% w Whole number
% ---------------------------------------------------------------------
w         = cfg_const;
w.tag     = 'w';
w.name    = 'Whole number';
w.val = {true};
% ---------------------------------------------------------------------
% i Integer
% ---------------------------------------------------------------------
i         = cfg_const;
i.tag     = 'i';
i.name    = 'Integer';
i.val = {true};
% ---------------------------------------------------------------------
% r Real number
% ---------------------------------------------------------------------
r         = cfg_const;
r.tag     = 'r';
r.name    = 'Real number';
r.val = {true};
% ---------------------------------------------------------------------
% strtype Type of output variable
% ---------------------------------------------------------------------
strtype         = cfg_choice;
strtype.tag     = 'strtype';
strtype.name    = 'Type of output variable';
strtype.values  = {s e n w i r };
% ---------------------------------------------------------------------
% any Any file
% ---------------------------------------------------------------------
any         = cfg_const;
any.tag     = 'any';
any.name    = 'Any file';
any.val = {true};
% ---------------------------------------------------------------------
% batch Batch file
% ---------------------------------------------------------------------
batch         = cfg_const;
batch.tag     = 'batch';
batch.name    = 'Batch file';
batch.val = {true};
% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_const;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.val = {true};
% ---------------------------------------------------------------------
% image Image(s)
% ---------------------------------------------------------------------
image         = cfg_const;
image.tag     = 'image';
image.name    = 'Image(s)';
image.val = {true};
% ---------------------------------------------------------------------
% mat MATLAB .mat file
% ---------------------------------------------------------------------
mat         = cfg_const;
mat.tag     = 'mat';
mat.name    = 'MATLAB .mat file';
mat.val = {true};
% ---------------------------------------------------------------------
% mesh Mesh file
% ---------------------------------------------------------------------
mesh         = cfg_const;
mesh.tag     = 'mesh';
mesh.name    = 'Mesh file';
mesh.val = {true};
% ---------------------------------------------------------------------
% nifti NIfTI file
% ---------------------------------------------------------------------
nifti         = cfg_const;
nifti.tag     = 'nifti';
nifti.name    = 'NIfTI file';
nifti.val = {true};
% ---------------------------------------------------------------------
% xml XML File
% ---------------------------------------------------------------------
xml         = cfg_const;
xml.tag     = 'xml';
xml.name    = 'XML File';
xml.val = {true};
% ---------------------------------------------------------------------
% filter Type of output file
% ---------------------------------------------------------------------
filter         = cfg_choice;
filter.tag     = 'filter';
filter.name    = 'Type of output file';
filter.values  = {any batch dir image mat mesh nifti xml };
% ---------------------------------------------------------------------
% outputs Outputs
% ---------------------------------------------------------------------
outputs         = cfg_repeat;
outputs.tag     = 'outputs';
outputs.name    = 'Outputs';
outputs.values  = {strtype filter };
outputs.num     = [0 Inf];
outputs.forcestruct = true;
% ---------------------------------------------------------------------
% fun Function to be called
% ---------------------------------------------------------------------
fun         = cfg_entry;
fun.tag     = 'fun';
fun.name    = 'Function to be called';
fun.strtype = 'f';
fun.num     = [1  1];
% ---------------------------------------------------------------------
% call_matlab Call MATLAB function
% ---------------------------------------------------------------------
call_matlab         = cfg_exbranch;
call_matlab.tag     = 'call_matlab';
call_matlab.name    = 'Call MATLAB function';
call_matlab.val     = {inputs outputs fun };
call_matlab.prog = @(job)cfg_run_call_matlab('run',job);
call_matlab.vout = @(job)cfg_run_call_matlab('vout',job);
% ---------------------------------------------------------------------
% cfg_basicio BasicIO
% ---------------------------------------------------------------------
cfg_basicio         = cfg_choice;
cfg_basicio.tag     = 'cfg_basicio';
cfg_basicio.name    = 'BasicIO';
cfg_basicio.help    = {'This toolbox contains basic input and output functions. The "Named Input" functions can be used to enter values or file names. These inputs can then be passed on to multiple modules, thereby ensuring all of them use the same input value. Some basic file manipulation is implemented in "Change Directory", "Make Directory", "Move Files". Lists of files can be filtered or splitted into parts using "File Set Filter" and "File Set Split". Output values from other modules can be written out to disk or assigned to MATLAB workspace.'};
cfg_basicio.values  = {file_move cfg_cd cfg_mkdir cfg_named_dir cfg_named_file file_fplist file_filter cfg_file_split cfg_named_input load_vars cfg_save_vars subsrefvar cfg_assignin runjobs call_matlab };
