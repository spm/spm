function item = cfg_const(varargin)

% This is the const configuration item class
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - 1x1 cell array
%    * check - (optional) function handle to implement configuration
%              specific subsasgn checks. This is not used here, since a
%              const should not change.
%    * help  - help text
% GUI/job manager fields
%    * expanded
%    * hidden
% All fields are inherited from the generic configuration item class.
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest
%    * all_set
%
% * 'const' - Constant value
%   - required fields: 'type', 'name', 'tag', 'val'
%   - optional fields: 'help'
%
%   The resulting data structure simply contains the contents of
%   val{1}.
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
% $Id: cfg_const.m 4166 2011-01-17 15:06:41Z volkmar $

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
    case {1,2,3,4,5}
        gitem = cfg_item(varargin{1:pnargin});
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
item = class(struct([]), myclass, gitem, cfg_leaf);
if nargin > mxpnargin
    item.cfg_item.val = varargin{mxpnargin+1};
    mxpnargin = mxpnargin+1;
end;
