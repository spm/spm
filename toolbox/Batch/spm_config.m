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
% * 'menu'  - Data entry by menu
%   - required fields: 'type', 'name', 'tag', 'labels', 'values'
%   - optional fields: 'val', 'def', 'help'
%
% * 'entry' - Entry by typing
%   - required fields: 'type', 'name', 'tag', 'strtype', 'num'
%   - optional fields: 'val', 'def', 'extras', 'help'
%
% * 'files' - Entry by file selection
%   - required fields: 'type', 'name', 'tag', 'filter', 'num'
%   - optional fields: 'val', 'def', 'help'
%
% * 'branch' - Branch of the tree structure
%   - required fields: 'type', 'name', 'tag', 'val'
%   - optional fields: 'prog', 'vfiles', 'check', 'modality', 'help'
%
% * 'choice' - A choice of ways of changing the tree structure
%   - required fields: 'type', 'name', 'tag', 'values'
%   - optional fields: 'val', 'help'
%
% * 'repeat' - Repeated kids in the tree structure
%   - required fields: 'type', 'name', 'tag'^, 'values'
%   - optional fields: 'help'
%
%     ^ except when length(values)==1
%
% Meanings of the fields:
% 'type'
% Indicates the type of each node in the tree structure (see above)
%
% 'name'
% A string containing a user readable name for the node.  This can
% contain spaces etc.
%
% 'tag'
% A name for the node that is used in the data structure saved by
% the job-system. It is not used for nodes of type 'repeat' with
% only one possible 'value'.
%
% 'help'
% Contains information for the user about what the node represents.
% It can be a string, or a cell array of strings. spm_justify.m
% can be used for aesthetic reasons.
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
% A filter for file selection - e.g. '*.mat'
%
% 'strtype'
% The type of values that are entered by typing.  e.g. 'e' for evaluated.
% The valid value types are:
%   's'   string
%   's+'  multi-line string
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
% to elements of a 'dim' field.
%
% 'prog'
%
% 'vfiles'
%
% 'check'
%
% 'modality'
% 
% 'dim' - obsolete.
% Dimensions.  This one is tricky to explain, but it allows dimensions
% to be shared between values that are descendents in the tree
% structure.  An element of a 'num' of e.g. -3, refers to dim(3), where
% dim is accumulated from all the dim fields in the ancestry of the
% node containing the num of -3.
%_______________________________________________________________________
%
% The easiest way to figure this stuff out is to look at examples in the
% spm_config_*.m files.
%_______________________________________________________________________
% %W% %E%

%persistent cnf
%if ~isempty(cnf),
%	vals = cnf;
%	return;
%end;

w       = spm_jobman('HelpWidth');

tempo   = struct('type','repeat','name','Temporal','tag','temporal',...
          'values',{{spm_config_slice_timing}},...
          'help',{spm_justify(w,'Temporal pre-processing functions.')});

spat    = struct('type','repeat','name','Spatial','tag','spatial',...
          'values',{{spm_config_realign,spm_config_realign_and_unwarp,spm_config_coreg,spm_config_segment,spm_config_norm,spm_config_smooth}},...
          'help',{spm_justify(w,'Various spatial and other pre-processing functions.')});
stat    = struct('type','repeat','name','Stats','tag','stats',...
          'values',{{spm_config_fmri_stats,spm_config_contrasts}},...
          'help',{spm_justify(w,'Various analysis utilities.')});
utils   = struct('type','repeat','name','Util','tag','util',...
          'values',{{spm_config_display,spm_config_checkreg,spm_config_imcalc,...
                     spm_config_dicom,spm_config_minc,spm_config_ecat}},...
          'help',{spm_justify(w,'Various useful tools.')});
ui      = spm_config_ui;

vals.type   = 'repeat';
vals.tag    = 'jobs';
vals.name   = 'SPM Jobs';
vals.values = {tempo,spat,stat,utils,ui};

[f,unused]  = spm_list_files(fullfile(spm('Dir'),'toolbox'),'*_config_*.m');
if ~isempty(f),
	tools      = struct('type','repeat','tag','tools','name','Tools','values',{{}});
	tools.help = spm_justify(w,'Other tools');
	for i=1:size(f,1),
		[unused,fl,unused] = fileparts(deblank(f(i,:)));
		if exist(fl,'file')~=2
			path(path,fullfile(spm('Dir'),'toolbox'));
		end;
		tools.values{i} = feval(fl);
	end;
	vals.values = {vals.values{:},tools};
end;


p1 = spm_justify(w,...
'The current list of jobs, which is represented as a tree-structure.',...
'Double-clicking can expand/contract items of the tree (marked with +/-)',...
'for visualisation.',...
'Items marked with X still require some values to be set before the job can',...
'be run, although an incompletely specified job can still be saved and loaded.');

p2 = spm_justify(w,...
'These are the options avaliable for the currently highlighted item.',...
'Changing the list of jobs is done by clicking on an option in the menu.',...
'Items can created, replicated or removed, allowing',...
'the processing stream to be modified.',...
'Values are also modified or entered via this panel.',...
'This is either by specifying values as text, selecting a menu option,',...
'or by file selection.');

p3 = spm_justify(w,...
'This panel shows the current value of the highlighted item (where relevant).');

p4 = spm_justify(w,...
'Jobs can be saved and loaded at a later time, either as',...
'XML or Matlab .mat files.  The format depends on the extension you give the filename.',...
'XML files can be loaded into Matlab via "loadxml", modified, and',...
'saved again by "savexml", whereas "load" and "save" can be used for Matlab .mat',...
'files. Incomplete jobs can be loaded or saved, but the specification needs to be',...
'complete for a job to be run.');

p5 = spm_justify(w,...
'This panel provides information about the meaning of the current item.');

vals.help   = {...
'* Top Left Panel',p1{:},'',...
'* Top Right Panel',p2{:},'',...
'* Centre Right Panel',p3{:},'',...
'* Save, Load & Run',p4{:},'',...
'* Bottom Panel',p5{:},''...
};

cnf = vals;
