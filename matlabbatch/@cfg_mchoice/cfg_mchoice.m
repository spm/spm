function item = cfg_mchoice(varargin)

% This is the choice configuration item class
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - 1x1 cell array of cfg_items (not set initially)
%    * check - (optional) function handle to implement configuration
%              specific subsasgn checks based on the harvested subtree
%              rooted at this node
%    * help  - help text
% GUI/job manager fields
%    * expanded
%    * hidden
% All fields are inherited from the generic configuration item class.
% Added fields
%    * values
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest     - a struct with a single field (see below)
%    * all_set     - returns all_set(item.val)
%
% Output in Job Structure (harvest)
% =================================
% The resulting data structure is a struct with a single field.  The
% name of the field is given by the 'tag' of the specified value.
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
% $Id: cfg_mchoice.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

myclass = mfilename;
% Get local fields and defaults from private/mysubs_fields
[fn, defs] = mysubs_fields;

if nargin == 1
    if isstruct(varargin{1})
        % assume input is a struct to be converted back into a class
        % return with error if this does not work
        if numel(fieldnames(varargin{1})) == numel(fn)+2 && ...
                all(isfield(varargin{1}, [fn(:)', {'cfg_item' 'cfg_intree'}]))
            gitem = varargin{1}.cfg_item;
            sitem = rmfield(varargin{1},{'cfg_item','cfg_intree'});
            item  = class(sitem, myclass, gitem, cfg_intree);
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
    case {1,2,3,4,5}
        gitem = cfg_item(varargin{1:pnargin});
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
for k=1:numel(fn)
    sitem.(fn{k})=defs{k};
end;
item = class(sitem, myclass, gitem, cfg_intree);
% set additional fields (if any) - field order as in mysubs_fields
for k = 1:min(numel(fn),nargin-mxpnargin)
    item.(fn{k}) = varargin{k+mxpnargin};
end;
