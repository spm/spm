function dep = cfg_dep(varargin)

% This is the configuration dependency class
%
% Data structure
% ==============
% Description fields
%    * sname        - display name of dependency source
%    * src_exbranch - subsref/subsasgn struct referencing the dependency
%                     source exbranch
%    * src_output   - subsref/subsasgn struct referencing the dependency
%                     source output item
%    * tname        - display name of dependency target
%    * tgt_exbranch - subsref/subsasgn struct referencing the dependency
%                     target exbranch in the config tree
%    * tgt_input    - subsref/subsasgn struct referencing the dependency
%                     target item in the config tree
%    * tgt_spec     - an cfg_findspec that can be used to restrict matches
%                     of this dependency to certain input items - this will
%                     be checked in subsasgn checks for cfg_items. Defaults
%                     to {} - match all cfg_items.
%    * jtsubs       - subsref/subsasgn struct referencing the dependency
%                     target item in the job tree (this is currently not
%                     used and may be removed in future)
%
% Public Methods
% ==============
%
% Public internal Methods
% =======================
%    * subsasgn
%    * subsref
%    * display
%    * disp
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_dep.m 4863 2012-08-27 08:09:23Z volkmar $

rev = '$Rev: 4863 $'; %#ok

dep = class(struct('tname','Target', ...
                   'tgt_exbranch', struct('type',{},'subs',{}), ...
                   'tgt_input', struct('type',{},'subs',{}), ...
                   'tgt_spec', {{}}, ...
                   'jtsubs', struct('type',{},'subs',{}), ...
                   'sname','Source', ...
                   'src_exbranch', struct('type',{},'subs',{}), ...
                   'src_output', struct('type',{},'subs',{})), ...
            'cfg_dep');
switch nargin
    case 0
        return;
    case 1
        if isa(varargin{1},'cfg_dep')
            dep = varargin{1};
        else
            dep.sname                          = varargin{1}; 
        end;
    case 2
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
    case 3
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
        dep.src_output(1:numel(varargin{3}))   = varargin{3};
    case 4
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
        dep.src_output(1:numel(varargin{3}))   = varargin{3};
        dep.tgt_spec                           = varargin{4};
    case 5
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
        dep.src_output(1:numel(varargin{3}))   = varargin{3};
        dep.tgt_spec                           = varargin{4};
        dep.jtsubs(1:numel(varargin{5}))       = varargin{5};
    case 6
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
        dep.src_output(1:numel(varargin{3}))   = varargin{3};
        dep.tgt_spec                           = varargin{4};
        dep.jtsubs(1:numel(varargin{5}))       = varargin{5};
        dep.tname                              = varargin{6};
    case 7
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
        dep.src_output(1:numel(varargin{3}))   = varargin{3};
        dep.tgt_spec                           = varargin{4};
        dep.jtsubs(1:numel(varargin{5}))       = varargin{5};
        dep.tname                              = varargin{6};
        dep.tgt_exbranch(1:numel(varargin{7})) = varargin{7};
    case 8
        dep.sname                              = varargin{1};
        dep.src_exbranch(1:numel(varargin{2})) = varargin{2};
        dep.src_output(1:numel(varargin{3}))   = varargin{3};
        dep.tgt_spec                           = varargin{4};
        dep.jtsubs(1:numel(varargin{5}))       = varargin{5};
        dep.tname                              = varargin{6};
        dep.tgt_exbranch(1:numel(varargin{7})) = varargin{7};
        dep.tgt_input(1:numel(varargin{8}))    = varargin{8};
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
