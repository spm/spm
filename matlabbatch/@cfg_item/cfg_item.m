function item = cfg_item(varargin)

% This is the generic configuration item class, from which all other
% classes are derived. 
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - (optional) val field: cell array
%    * check - (optional) function handle to implement configuration
%              specific checks based on the harvested subtree rooted at
%              this node. It will be evaluated during harvest if all
%              dependencies in the harvested subtree are resolved and all
%              val's are set. 
%              This function should return an empty string on success and
%              a string explaining why it failed otherwise. 
%    * help  - help text
% GUI/job manager fields
%    * id   
%    * expanded
%    * hidden
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%                    No validity check performed here, this needs to be
%                    added in child class method.
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest     - returns item.val{1}, or '<UNDEFINED>' if empty
%    * all_set     - returns ~isempty(item.val)
%
% Public internal Methods
% =======================
%    * subsasgn
%    * subsref
%    * display
%    * disp
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
% $Id: cfg_item.m 1195 2008-03-07 21:51:49Z volkmar $

rev = '$Rev: 1195 $';

myclass = mfilename;
% Get local fields and defaults from private/mysubs_fields
[fn defs] = mysubs_fields;

if nargin == 1
    if isstruct(varargin{1})
        % assume input is a struct to be converted back into a class
        % return with error if this does not work
        if numel(fieldnames(varargin{1})) == numel(fn) && all(isfield(varargin{1}, fn))
            item = class(varargin{1}, myclass);
            return;
        else
            error('matlabbatch:reclassify', ['Don''t know how to convert this ' ...
                            'into class ''%s''.'], myclass);
        end;
    end;
    if isa(varargin{1}, myclass)
        item = varargin{1};
        return;
    end;
end;

item = class(struct('name','Config Item', ...
    'tag','generic', 'val',{{}}, 'check',{''}, 'help',{{''}}, 'id',[], 'expanded',0, 'hidden',0),...
    'cfg_item');
switch nargin
    case 0
        return;
    case 1
        item.name = varargin{1};
    case 2
        item.name = varargin{1};
        item.tag  = varargin{2};
    case 3
        item.name  = varargin{1};
        item.tag   = varargin{2};
        item.check = varargin{3};
    case 3
        item.name  = varargin{1};
        item.tag   = varargin{2};
        item.check = varargin{3};
        item.help  = varargin{4};
    otherwise
        error('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
