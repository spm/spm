function vals = spm_config
% Configuration file for SPM jobs
%
% A configuration structure consists of a tree of data-structures.
% Each node is a structure, which must contain a field called "type".
% Allowable types are:
% * 'const' - Constant value
%   - required fields: 'type', 'name', 'tag', 'val'
%   - optional fields: 'help'
%
%   The resulting data structure simply contains the contents of
%   val{1}.
%
%
% * 'menu'  - Data entry by menu
%   - required fields: 'type', 'name', 'tag', 'labels', 'values'
%   - optional fields: 'val', 'def', 'help'
%
%   The resulting data structure simply contains the contents of
%   val{1}, which corresponds to the element of values selected
%   by the user.
%
%
% * 'entry' - Entry by typing
%   - required fields: 'type', 'name', 'tag', 'strtype', 'num'
%   - optional fields: 'val', 'def', 'extras', 'help'
%
%   The resulting data structure simply contains the contents of
%   val{1}, which is what the user specified.
%
%
% * 'files' - Entry by file selection
%   - required fields: 'type', 'name', 'tag', 'filter', 'num'
%   - optional fields: 'val', 'def', 'help'
%
%   The resulting data structure simply contains a cell array
%   of filenames (from val{1} ).
%
%
% * 'branch' - Branch of the tree structure
%   - required fields: 'type', 'name', 'tag', 'val'
%   - optional fields: 'prog', 'vfiles', 'check', 'modality', 'help'
%
%   The resulting data structure is a struct, with fieldnames according
%   to the 'tag's of the child nodes.
%
%
% * 'choice' - A choice of ways of changing the tree structure
%   - required fields: 'type', 'name', 'tag', 'values'
%   - optional fields: 'val', 'help'
%
%   The resulting data structure is a struct with a single field.  The
%   name of the field is given by the 'tag' of the specified value.
%
%
% * 'repeat' - Repeated kids in the tree structure
%   - required fields: 'type', 'name', 'tag', 'values'
%   - optional fields: 'num', 'help'
%
%   If the number of elements in the 'values' field is greater than
%   one, then the resulting data structure is a cell array.  Each
%   element of the cell array is a struct with a single field, where
%   the name of the field is given by the 'tag' of the child node.
%
%   If the 'values' field only contains one element, which is a 'branch',
%   then the data structure is a struct array (in which case the 'tag'
%   of the current node can be ignored).
%
%   If the 'values' field only contains one element, which is not a branch,
%   then the data structure is a cell array, where each element is the value
%   of the child node ((in which case the 'tag' of the current node can be
%   ignored).
%
%   If the 'num' field is present, it is a 2-vector [min max] which
%   limits the occurence of the repeated substructure(s). This can
%   be used to check inputs whether there is at least one repeated
%   entry.
%
%
%
% Meanings of the fields:
% 'type'
% Indicates the type of each node in the tree structure (see above)
%
% 'name'
% A string containing a user readable name for the node.  This can
% contain spaces etc. This is what is used to prompt the user for
% entries.
%
% 'tag'
% A name for the node that is used in the data structure saved by
% the job-system. It is not used for nodes of type 'repeat' with
% only one possible 'value'.
%
% 'help'
% Contains information for the user about what the node represents.
% If it is a cell array, then each element will be treated as a
% paragraph, and justified accordingly. See spm_justify.
%
% 'values'
% A cell vector of possible values.  It has different meanings in
% different contexts.  For 'menu' types, it is a list of values that
% a variable can have.  For 'choice' and 'repeat' types, each
% element is a node.
%
% 'val'
% The cell vector containing the current value(s) of that node.  For
% 'branch', 'choice' and 'repeat' types, the value(s) are nodes. For
% 'const', 'menu', 'entry' and 'files' types, they contain a value
% for that variable.
%
% 'def'
% A string containing the name of an element in the 'defaults'
% structure (see spm_defaults.m). Undefined 'val's are assigned
% the value of this default.
%
% 'labels'
% A cell array of strings containing names for the different values.
% It is used by nodes of type 'menu'.
%
% 'filter'
% A filter for file selection.  See spm_select for filter usage.
%
% 'strtype'
% The type of values that are entered by typing.  e.g. 'e' for evaluated.
% The valid value types are:
%   's'   string
%   'e'   evaluated expression
%   'n'   natural number (1..n)
%   'w'   whole number (0..n)
%   'i'   integer
%   'r'   real number
%   'c'   indicator vector (e.g., 0101... or abab...)
%   'x'   contrast matrix
%   'p'   permutation
%
% 'extras'
% Extra arguments for obtaining values that are typed in (e.g. min and
% max values).
%
% 'num'
% A one or two element vector containing the size of the required
% value for the node.  For 'entry' nodes, it is the size of an
% array.  For 'files' nodes, it is the min and max number of files.
% Negative values have special meaning in both cases.  These refer
% to elements of a 'dim' field.  See spm_select for more.
%
% 'prog'
% This is a function handle that may be invoked.  The tree structure
% of user specified things (starting at the point with this field) is
% collected together and passed as the first (and only) argument to this
% function via "feval". It is executed when the data structure is run
% (usually when the "Run" button is hit).
% e.g. opts.prog = @some_function;
% 
% 'vfiles'
% Another function handle.  This is usually used to produce a list of
% files that would be created when the program is run.  This list is in the
% form of a cell array of strings.  When creating batch jobs, it is useful
% to be able to select files that have not yet been created (see spm_select).
% As in the case of the 'prog' field, the user-specified things are collected
% together and passed as the first argument.
% e.g. opts.prog = @some_function;
%
% 'check'
% Another function handle.  Sometimes it is advisable to do some checking of
% what the user has entered (to ensure dimensions match etc).  Note that
% it is not a good idea for the checking to have to open files etc, as some
% of the files that have been selected may not yet exist.  Again, the user
% specified data structure is collected together and passed as the first
% argument.  The output is either empty, or it is a string that explains
% what the problem is.
% e.g. opts.prog = @some_function;
%
% 'modality'
% A cell array of strings indicating which modalities this branch of the tree
% structure should be visible for.  If there is no modality field, then it is
% visible for all.
% e.g. opts.modality = {'PET','FMRI'};
%      opts.modality = {'EEG'};
%
%_______________________________________________________________________
%
% The easiest way to figure this stuff out is to look at working
% examples in the spm_config_*.m files.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config.m 232 2005-09-15 19:02:59Z john $


tempo   = struct('type','repeat','name','Temporal','tag','temporal',...
          'values',{{spm_config_slice_timing,spm_config_eeg_filter}},...
          'help',{'Temporal pre-processing functions.'});

spat    = struct('type','repeat','name','Spatial','tag','spatial',...
          'values',{{spm_config_realign,spm_config_realign_and_unwarp,spm_config_coreg,...
                     spm_config_preproc,spm_config_norm,spm_config_smooth}},...
          'help',{'Various spatial and other pre-processing functions.'});
stat    = struct('type','repeat','name','Stats','tag','stats',...
          'values',{{spm_config_fmri_spec,spm_config_fmri_est,spm_config_contrasts}},...
          'help',{'Various analysis utilities.'});
utils   = struct('type','repeat','name','Util','tag','util',...
          'values',{{spm_config_display,spm_config_checkreg,spm_config_imcalc,...
                     spm_config_dicom,spm_config_minc,spm_config_ecat,...
                     spm_config_runbatch,...
                     spm_config_cd,spm_config_mkdir,spm_config_defs,spm_config_ui}},...
          'help',{'Various useful tools.'});

vals.type   = 'repeat';
vals.tag    = 'jobs';
vals.name   = 'SPM Jobs';
vals.values = {tempo,spat,stat,utils};

[f,unused]  = spm_list_files(fullfile(spm('Dir'),'toolbox'),'*_config_*.m');
if ~isempty(f),
	tools      = struct('type','repeat','tag','tools','name','Tools','values',{{}});
	tools.help = {'Other tools'};
	for i=1:size(f,1),
		[unused,fl,unused] = fileparts(deblank(f(i,:)));
		if exist(fl,'file')~=2
			path(path,fullfile(spm('Dir'),'toolbox'));
		end;
		tools.values{i} = feval(fl);
	end;
	vals.values = {vals.values{:},tools};
end;


p1 = [...
'The current list of jobs, which is represented as a tree-structure. ',...
'Double-clicking can expand/contract items of the tree (marked with +/-) ',...
'for visualisation. ',...
'Items marked with X still require some values to be set before the ',...
'job can be run, although an incompletely specified job can still be ',...
'saved and loaded.'];

p2 = [...
'These are the options avaliable for the currently highlighted item. ',...
'Changing the list of jobs is done by clicking on an option in the menu. ',...
'Items can created, replicated or removed, allowing ',...
'the processing stream to be modified. ',...
'Values are also modified or entered via this panel. ',...
'This is either by specifying values as text, selecting a menu option, ',...
'or by file selection.'];

p3 = ...
'This panel shows the current value of the highlighted item (where relevant).';

p4 = [...
'Jobs can be saved and loaded at a later time, either as ',...
'XML or Matlab .mat files.  The format depends on the extension you ',...
'give the filename. XML files can be loaded into Matlab via "loadxml", ',...
'modified, and saved again by "savexml", whereas "load" and "save" ',...
'can be used for Matlab .mat files. Incomplete jobs can be loaded or ',...
'saved, but the specification needs to be ',...
'complete for a job to be run.'];

p5 = ...
'This panel provides information about the meaning of the current item.';

fg = [...
'/*\begin{figure} ',...
'\begin{center} ',...
'\epsfig{file=ui1,width=70mm} \epsfig{file=ui2,width=70mm} ',...
'\epsfig{file=ui3,width=70mm} \epsfig{file=ui4,width=70mm}',...
'\end{center} ',...
'\caption{The SPM5 user interface. ',...
'Top left: The usual user-interface.  Top right: The Defaults user-interface. ',...
'Bottom left: The file selector (click the (?) button for more information about ',...
'filtering filenames, or selecting individual volumes within a 4D file. ',...
'Bottom right: more online help can be obtained via the main help button.} ',...
'\end{figure} */'];
vals.help   = {...
'%* Top Left Panel','/*\subsection*{Top Left Panel}*/',p1,'',...
'%* Top Right Panel','/*\subsection*{Top Right Panel}*/',p2,'',...
'%* Centre Right Panel','/*\subsection*{Centre Right Panel}*/',p3,'',...
'%* Save, Load & Run','/*\subsection*{Save, Load \& Run}*/',p4,'',...
'%* Bottom Panel','/*\subsection*{Bottom Panel}*/',p5,fg...
};

