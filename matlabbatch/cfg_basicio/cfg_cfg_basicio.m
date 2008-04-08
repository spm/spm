function cfg_basicio = cfg_cfg_basicio
% 'BasicIO' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2008-03-12 09:23:12.
% ---------------------------------------------------------------------
% files Files to move/delete
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'Files to move/delete';
files.help    = {'These files will be moved/deleted.'};
files.filter = 'any';
files.ufilter = '.*';
files.num     = [1 Inf];
% ---------------------------------------------------------------------
% target Target Directory
% ---------------------------------------------------------------------
target         = cfg_files;
target.tag     = 'target';
target.name    = 'Target Directory';
target.help    = {'Files will be moved to the specified directory. If no directory is specified, then files will be deleted instead.'};
target.filter = 'dir';
target.ufilter = '.*';
target.num     = [0 1];
% ---------------------------------------------------------------------
% file_move Move/Delete Files
% ---------------------------------------------------------------------
file_move         = cfg_exbranch;
file_move.tag     = 'file_move';
file_move.name    = 'Move/Delete Files';
file_move.val     = {files target };
file_move.help    = {'Move or delete files.'};
file_move.prog = @cfg_run_file_move;
file_move.vout = @cfg_vout_move_file;
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
name.num     = [1 Inf];
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
name.num     = [1 Inf];
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
name.num     = [1 Inf];
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
filter.help    = {'A regular expression to filter files (applied after filtering for ''Typ'').'};
filter.strtype = 's';
filter.num     = [1 Inf];
% ---------------------------------------------------------------------
% file_fplist File Selector (Batch Mode)
% ---------------------------------------------------------------------
file_fplist         = cfg_exbranch;
file_fplist.tag     = 'file_fplist';
file_fplist.name    = 'File Selector (Batch Mode)';
file_fplist.val     = {dir filter };
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
typ.num     = [1 Inf];
% ---------------------------------------------------------------------
% filter Filter
% ---------------------------------------------------------------------
filter         = cfg_entry;
filter.tag     = 'filter';
filter.name    = 'Filter';
filter.help    = {'A regular expression to filter files (applied after filtering for ''Typ'').'};
filter.strtype = 's';
filter.num     = [1 Inf];
% ---------------------------------------------------------------------
% frames Frames
% ---------------------------------------------------------------------
frames         = cfg_entry;
frames.tag     = 'frames';
frames.name    = 'Frames';
frames.help    = {'Number of frames to address in 4D NIfTI images. Ignored if empty.'};
frames.strtype = 's';
frames.num     = [0 Inf];
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
name.num     = [1 Inf];
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
index.num     = [1 Inf];
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
name.num     = [1 Inf];
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
% name Output Filename
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Output Filename';
name.help    = {'Output filename without any directory names.'};
name.strtype = 's';
name.num     = [1 Inf];
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
vname.num     = [1 Inf];
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
saveasstruct.values{1} = logical(false);
saveasstruct.values{2} = logical(true);
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
% name Output Variable Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Output Variable Name';
name.check   = @cfg_check_assignin;
name.help    = {'Enter a valid MATLAB variable name. The contents of the input to "Output Item" will be assigned to a variable of this name in MATLAB workspace.'};
name.strtype = 's';
name.num     = [1 Inf];
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
                        'The value entered into "Output Item" will be assigned to a MATLAB workspace variable whose name is specified in "Output Variable Name". This can be useful to assess the results of computations with other MATLAB routines or for debugging.'
}';
cfg_assignin.prog = @cfg_run_assignin;
% ---------------------------------------------------------------------
% cfg_basicio BasicIO
% ---------------------------------------------------------------------
cfg_basicio         = cfg_repeat;
cfg_basicio.tag     = 'cfg_basicio';
cfg_basicio.name    = 'BasicIO';
cfg_basicio.help    = {'This toolbox contains basic input and output functions. The "Named Input" functions can be used to enter values or file names. These inputs can then be passed on to multiple modules, thereby ensuring all of them use the same input value. Some basic file manipulation is implemented in "Change Directory", "Make Directory", "Move Files". Lists of files can be filtered or splitted into parts using "File Set Filter" and "File Set Split". Output values from other modules can be written out to disk or assigned to MATLAB workspace.'};
cfg_basicio.values  = {file_move cfg_cd cfg_mkdir cfg_named_dir cfg_named_file file_fplist file_filter cfg_file_split cfg_named_input cfg_save_vars cfg_assignin };
cfg_basicio.num     = [0 Inf];
cfg_basicio.forcestruct = true;
