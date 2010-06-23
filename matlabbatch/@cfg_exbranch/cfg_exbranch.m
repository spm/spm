function item = cfg_exbranch(varargin)

% This is the exbranch configuration item class
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - 1xn cell array of cfg_items
%    * check - (optional) function handle to implement configuration
%              specific subsasgn checks based on the harvested subtree
%              rooted at this node
%    * help  - help text
% GUI/job manager fields
%    * expanded
%    * hidden
% All fields above are inherited from the branch configuration item class.
%    * prog
%    * vfiles
%    * modality
%    * vout  - function handle that generates sout struct
%    * sout  - source dependency description
%    * jout  - saved output (will be referenced by harvest of a dependency
%              target for dependency resolution at job runtime)
%    * tdeps - list where this branch is target of a dependency
%    * sdeps - list where this branch is source of a dependency
%    * chk   - field to save check status from cfg_item.check callbacks
%    * id    - id of this cfg_exbranch. This is used to reference the
%              cfg_exbranch in cfg_dep objects.
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest
%    * all_set
%
% * 'executable branch'  - See branch for details on inherited fields
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
% $Id: cfg_exbranch.m 3944 2010-06-23 08:53:40Z volkmar $

rev = '$Rev: 3944 $'; %#ok

myclass = mfilename;
% Get local fields and defaults from private/mysubs_fields
[fn defs] = mysubs_fields;

if nargin == 1
    if isstruct(varargin{1})
        % assume input is a struct to be converted back into a class
        % return with error if this does not work
        if numel(fieldnames(varargin{1})) == numel(fn)+1 && ...
                all(isfield(varargin{1}, [fn(:)' {'cfg_branch'}]))
            bitem = varargin{1}.cfg_branch;
            sitem = rmfield(varargin{1},{'cfg_branch'});
            item  = class(sitem, myclass, bitem);
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

mxpnargin = 5; % Max 5 arguments to parent initialisation
pnargin = min([nargin, mxpnargin]); 
switch nargin
    case 0
        bitem = cfg_branch;
    case 1
        if isa(varargin{1},myclass)
            item = varargin{1};
            return;
        elseif isstruct(varargin{1})
            % assume input is a struct to be converted back into a class
            % no fieldname checking performed here
            bitem = varargin{1}.cfg_branch;
            sitem = rmfield(varargin{1},'cfg_branch');
            item  = class(sitem, myclass, bitem);
            return;
        else
            bitem = cfg_branch(varargin{1});
        end;
    case {2,3,4,5,6,7,8,9,10}
        bitem = cfg_branch(varargin{1:pnargin});
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
for k=1:numel(fn)
    sitem.(fn{k})=defs{k};
end;
item = class(sitem, myclass, bitem);
% set additional fields (if any) - field order as in mysubs_fields
for k = 1:min(numel(fn),nargin-mxpnargin)
    item.(fn{k}) = varargin{k+mxpnargin};
end;
