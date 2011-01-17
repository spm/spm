function item = cfg_files(varargin)

% This is the file configuration item class
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - 1x1 cell array
%    * check - (optional) function handle to implement configuration
%              specific subsasgn checks based on the harvested subtree
%              rooted at this node
%    * help  - help text
% GUI/job manager fields
%    * expanded
%    * hidden
% All fields above are inherited from the generic configuration item class.
%    * filter  - default '.*'
%    * num     - default [0 Inf]
%    * dir     - default ''
%    * ufilter - default '.*'
%    * def
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest     - returns item.val{1}, or '<UNDEFINED>' if empty, see below
%    * all_set     - returns ~isempty(item.val), checks numel(item.val{1})
%                    against item.num
%
% Subscript Assignment Checks
% ===========================
% Values assigned to the .num field must be a 2-vector.
% Values assigned to the .val{1} field must be either
% - empty
% - an array of cfg_dep objects
% - a cell string of file names.
% In the latter case, the cell string will be filtered using 
% cfg_getfile('filter', item.val{1}, item.filter, '.*', 1:inf) and only
% files matching item.filter will be assigned.
%
% GUI Input
% =========
% The GUI uses 
% cfg_getfile(item.num, item.filter, item.name, item.val{1}, '.', item.ufilter)
% to select files. The filter in item.filter can not be overridden by the
% GUI.
%
% Output in Job Structure (harvest)
% =================================
% cfg_files uses cfg_item/harvest. If multiple dependencies are present
% and all can be resolved, the result will be a cell string containing a
% concatenated list of files.
%
% The layout of the configuration tree and the types of configuration items
% have been kept compatible to a configuration system and job manager
% implementation in SPM5 (Statistical Parametric Mapping, Copyright (C)
% 2005 Wellcome Department of Imaging Neuroscience). This code has been
% completely rewritten based on an object oriented model of the
% configuration tree.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_files.m 4166 2011-01-17 15:06:41Z volkmar $

rev = '$Rev: 4166 $'; %#ok

myclass = mfilename;
% Get local fields and defaults from private/mysubs_fields
[fn defs] = mysubs_fields;

if nargin == 1
    if isstruct(varargin{1})
        % assume input is a struct to be converted back into a class
        % return with error if this does not work
        if numel(fieldnames(varargin{1})) == numel(fn)+2 && ...
                all(isfield(varargin{1}, [fn(:)', {'cfg_item' 'cfg_leaf'}]))
            gitem = varargin{1}.cfg_item;
            sitem = rmfield(varargin{1},{'cfg_item', 'cfg_leaf'});
            item  = class(sitem, myclass, gitem, cfg_leaf);
            return;
        else
            cfg_message('matlabbatch:constructor:reclassify', ['Don''t know how to convert this ' ...
                            'into class ''%s''.'], myclass);
        end;
    end;
    if isa(varargin{1},myclass)
        item = varargin{1};
        return;
    end;
end;

mxpnargin = 4; % Max 4 arguments to parent initialisation
pnargin = min([nargin,mxpnargin]);
switch nargin
    case 0
        gitem = cfg_item;
    case {1,2,3,4,5,6,7,8,9,10}
        gitem = cfg_item(varargin{1:pnargin});
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
for k=1:numel(fn)
    sitem.(fn{k})=defs{k};
end;
item = class(sitem, myclass, gitem, cfg_leaf);
if nargin > mxpnargin
    item.cfg_item.val{1} = varargin{mxpnargin+1};
    mxpnargin = mxpnargin+1;
end;
% set additional fields (if any) - field order as in mysubs_fields
for k = 1:min(numel(fn),nargin-mxpnargin)
    item.(fn{k}) = varargin{k+mxpnargin};
end;
